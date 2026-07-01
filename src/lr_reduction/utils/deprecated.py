"""Decorator to mark functions as deprecated. Only needed for Python < 3.13.
For Python >= 3.13, use warnings.deprecated instead."""

import functools
import inspect
import warnings


def deprecated(message: str):
    """
    A decorator to mark functions as deprecated. Emits a warning when the function is invoked.
    """

    def decorator(func):
        @functools.wraps(func)
        def new_func(*args, **kwargs):
            warnings.simplefilter("always", DeprecationWarning)
            frame = inspect.currentframe()
            caller = frame.f_back if frame is not None else None
            module_name = caller.f_globals.get("__name__", func.__module__) if caller else func.__module__
            lineno = caller.f_lineno if caller else 1
            warnings.warn_explicit(
                message,
                category=DeprecationWarning,
                filename=module_name,
                lineno=lineno,
                module=module_name,
            )
            return func(*args, **kwargs)

        return new_func

    return decorator
