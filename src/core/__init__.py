"""
Core package

This package contains core application functionality:
- State management (AppState class and global instance)
- Logging utilities (Logger class and global instance)

"""

from .state import AppState, app_state
from .logger import (
    Logger, 
    logger, 
    # append_log, 
    # append_log_file, 
    # clear_logs
)

__all__ = [
    # State management
    'AppState',
    'app_state',
    
    # Logging
    'Logger',
    'logger',
    # 'append_log',
    # 'append_log_file',
    # 'clear_logs',
]