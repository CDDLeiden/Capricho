"""Module for rate limiting function calls through a decorator."""

import threading
import time
from functools import wraps

from ..logger import logger


def rate_limit(max_per_second=5):
    min_interval = 1.0 / max_per_second
    lock = threading.Lock()
    last_called = 0

    logger.debug(f"Initializing rate limiter: max {max_per_second} calls per second")

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            nonlocal last_called
            with lock:
                current_time = time.time()
                elapsed = current_time - last_called
                left_to_wait = min_interval - elapsed

                logger.trace(f"Function {func.__name__} called. Time since last call: {elapsed:.4f}s")

                if left_to_wait > 0:
                    logger.trace(f"Rate limit exceeded. Waiting for {left_to_wait:.4f}s")
                    time.sleep(left_to_wait)
                else:
                    logger.trace("Executing immediately")

                try:
                    logger.trace(f"Executing {func.__name__}")
                    ret = func(*args, **kwargs)
                    logger.trace(f"Finished executing {func.__name__}")
                except Exception as e:
                    logger.exception(f"Exception in {func.__name__}: {str(e)}")
                    raise

                last_called = time.time()
                logger.trace(f"Updated last_called to {last_called}")
            return ret

        return wrapper

    return decorator
