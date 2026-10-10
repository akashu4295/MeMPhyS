"""Start the MeMPhyS PySide6 desktop application."""

import sys

from src.config.constants import APP_FULL_NAME, LOG_FILE_PATH
from src.core import app_state, logger

try:
    from src.qt_app import main as qt_main
except ModuleNotFoundError as error:
    if error.name != "PySide6":
        raise
    qt_main = None


def main():
    if qt_main is None:
        print("PySide6 is not installed in the active Python environment.")
        print("Run with: conda run -n memphys_gui python gui.py")
        return 1

    log_handle = app_state.open_log_file(LOG_FILE_PATH, mode="a")
    if log_handle:
        logger.set_file_handle(log_handle)
    logger.info(f"Starting {APP_FULL_NAME}")
    try:
        return qt_main()
    except Exception as error:
        logger.log_exception(error, "Application startup failed")
        raise
    finally:
        logger.info("Application shutting down")
        app_state.cleanup()


if __name__ == "__main__":
    sys.exit(main())