"""Read-only source-data retrieval helpers for optional reference extractors.

Ordinary application runs never import this module or access the network.
ZIP members are decompressed as data; no deposited code is executed.
"""

from io import BytesIO, RawIOBase
import pickle
import re
from urllib.request import Request, urlopen


def fetch(url, limit=64*1024**2):
    with urlopen(Request(url, headers={"User-Agent": "qdjj_solver research applications"}), timeout=120) as response:
        data = response.read(limit+1)
    if len(data) > limit:
        raise ValueError("reference download exceeded the declared size limit")
    return data


def numeric_pickle(data):
    """Decode the legacy deposit's tuples of numeric NumPy arrays, with no arbitrary globals."""
    import numpy as np
    reconstruct = np.ndarray.__reduce__(np.array([]))[0]
    allowed = {("numpy", "ndarray"): np.ndarray, ("numpy", "dtype"): np.dtype,
               ("numpy.core.multiarray", "_reconstruct"): reconstruct,
               ("numpy._core.multiarray", "_reconstruct"): reconstruct,
               ("numpy.core.multiarray", "scalar"): np.float64(0).__reduce__()[0],
               ("numpy._core.multiarray", "scalar"): np.float64(0).__reduce__()[0]}

    class NumericUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            if (module, name) not in allowed:
                raise ValueError(f"unexpected global in numeric reference: {module}.{name}")
            return allowed[module, name]

    result = NumericUnpickler(BytesIO(data)).load()
    if not isinstance(result, (tuple, list)):
        raise ValueError("expected a tuple of numeric reference arrays")
    for item in result:
        array = np.asarray(item)
        if array.dtype.kind not in "fiu" or array.size > 1_000_000 or not np.all(np.isfinite(array)):
            raise ValueError("reference contains nonnumeric, nonfinite or oversized arrays")
    return result


class RemoteZip(RawIOBase):
    """Seekable HTTP Range view, suitable for zipfile.ZipFile.

    Verify Content-Range and the published ETag on every response. ZIP's CRC
    verifies each extracted member; this does not hash the entire remote ZIP.
    """

    def __init__(self, url, etag):
        self.url, self.etag = url, etag
        self.position = 0
        self.transferred = 0
        self.cache_start, self.cache = 0, b""
        self.size = None
        start, data = self._request("bytes=-65536")
        self.cache_start, self.cache = start, data

    def _request(self, byte_range):
        request = Request(self.url, headers={"Range": byte_range, "Accept-Encoding": "identity"})
        with urlopen(request, timeout=120) as response:
            match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
            if response.status != 206 or match is None:
                raise ValueError("reference server did not honor HTTP Range")
            if response.headers.get("ETag", "").strip('"') != self.etag:
                raise ValueError("reference archive ETag differs from the pinned publication")
            start, end, size = map(int, match.groups())
            if self.size is not None and size != self.size:
                raise ValueError("reference archive changed during retrieval")
            self.size = size
            if end-start+1 > 64*1024**2:
                raise ValueError("requested reference member is unexpectedly large")
            data = response.read(end-start+2)
            if len(data) != end-start+1:
                raise ValueError("incomplete HTTP Range response")
        self.transferred += len(data)
        return start, data

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        bases = {0: 0, 1: self.position, 2: self.size}
        if whence not in bases or bases[whence]+offset < 0:
            raise ValueError("invalid archive seek")
        self.position = bases[whence]+offset
        return self.position

    def read(self, size=-1):
        size = self.size-self.position if size < 0 else min(size, self.size-self.position)
        if size <= 0:
            return b""
        if not (self.cache_start <= self.position and self.position+size <= self.cache_start+len(self.cache)):
            end = min(self.size-1, self.position+max(size, 65536)-1)
            self.cache_start, self.cache = self._request(f"bytes={self.position}-{end}")
        start = self.position-self.cache_start
        self.position += size
        return self.cache[start:start+size]
