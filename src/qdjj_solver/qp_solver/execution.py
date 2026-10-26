"""Optional, cooperative execution hooks shared by fresh and prepared QP solves."""

from time import perf_counter


class CalculationCancelled(KeyboardInterrupt):
    """A requested cancellation, distinct from a numerical convergence failure."""


class ExecutionContext:
    """Cooperative QP progress and cancellation callbacks.

    ``progress(event)`` receives a dictionary with the stage, elapsed seconds and
    stage-specific details. Repeated events in one stage are separated by at least
    ``interval`` seconds. ``cancelled()`` is checked at every checkpoint and a true
    result raises ``CalculationCancelled``. Callback exceptions propagate to the caller.
    """

    def __init__(self, progress=None, cancelled=None, interval=0.1):
        self.progress = progress
        self.cancelled = cancelled
        self.interval = max(0., float(interval))
        self.started = perf_counter()
        self.last_event = -float("inf")
        self.last_stage = None

    def checkpoint(self, stage, **details):
        """Check cancellation and emit a progress event when the stage or interval changes."""
        # Cancellation is checked even when a progress event is throttled.
        if self.cancelled is not None and self.cancelled():
            raise CalculationCancelled("calculation cancelled")
        now = perf_counter()
        if self.progress is not None and (stage != self.last_stage or now-self.last_event >= self.interval):
            self.progress(dict(stage=stage, elapsed=now-self.started, **details))
            self.last_event, self.last_stage = now, stage
