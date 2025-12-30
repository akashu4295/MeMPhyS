"""
Solver monitoring module for MeMPhyS GUI

Monitors convergence data and updates plots in real-time
during solver execution.
"""

import os
import time
import threading
from typing import Optional
import numpy as np
import dearpygui.dearpygui as dpg

from src.config import CONVERGENCE_CSV, CONVERGENCE_UPDATE_INTERVAL
from src.core import logger, app_state
from src.utils import read_csv_file


class ConvergenceMonitor:
    """
    Monitors convergence CSV file and updates GUI plot in real-time
    
    Features:
    - Automatic file polling
    - Dynamic axis scaling
    - Thread-safe GUI updates
    - Graceful shutdown
    """
    
    def __init__(self, 
                 csv_file: str = CONVERGENCE_CSV,
                 series_tag: str = "conv_series",
                 x_axis_tag: str = "x_axis_conv",
                 y_axis_tag: str = "y_axis_conv",
                 update_interval: float = CONVERGENCE_UPDATE_INTERVAL):
        """
        Initialize the convergence monitor
        
        Args:
            csv_file: Path to convergence CSV file
            series_tag: DearPyGUI tag for line series
            x_axis_tag: DearPyGUI tag for x-axis
            y_axis_tag: DearPyGUI tag for y-axis
            update_interval: Update interval in seconds
        """
        self.csv_file = csv_file
        self.series_tag = series_tag
        self.x_axis_tag = x_axis_tag
        self.y_axis_tag = y_axis_tag
        self.update_interval = update_interval
        
        self.is_monitoring = False
        self.monitor_thread: Optional[threading.Thread] = None
        self.last_modified_time = 0
    
    def start(self):
        """Start monitoring in a background thread"""
        if self.is_monitoring:
            logger.warning("Convergence monitor is already running")
            return
        
        logger.info("Starting convergence monitor")
        self.is_monitoring = True
        app_state.convergence_monitor_running = True
        
        self.monitor_thread = threading.Thread(
            target=self._monitor_loop,
            daemon=True,
            name="ConvergenceMonitor"
        )
        self.monitor_thread.start()
    
    def stop(self):
        """Stop monitoring"""
        if not self.is_monitoring:
            return
        
        logger.info("Stopping convergence monitor")
        self.is_monitoring = False
        app_state.convergence_monitor_running = False
        
        # Wait for thread to finish (with timeout)
        if self.monitor_thread:
            self.monitor_thread.join(timeout=self.update_interval * 2)
        
        logger.info("Convergence monitor stopped")
    
    def _monitor_loop(self):
        """Main monitoring loop (runs in separate thread)"""
        logger.debug("Convergence monitor loop started")
        
        # Wait a bit for DearPyGUI to be fully running
        import time
        time.sleep(1)
        
        while self.is_monitoring:
            try:
                # Check if file exists and has been modified
                if self._should_update():
                    self._update_plot()
                
                # Sleep before next check
                time.sleep(self.update_interval)
            
            except Exception as e:
                # Only log if we're still supposed to be running
                if self.is_monitoring:
                    logger.debug(f"Error in convergence monitor: {e}")
                time.sleep(self.update_interval)
        
        logger.debug("Convergence monitor loop ended")
    
    def _should_update(self) -> bool:
        """
        Check if the convergence file should be read
        
        Returns:
            True if file exists and has been modified since last read
        """
        if not os.path.exists(self.csv_file):
            return False
        
        try:
            current_modified_time = os.path.getmtime(self.csv_file)
            
            if current_modified_time > self.last_modified_time:
                self.last_modified_time = current_modified_time
                return True
            
            return False
        
        except Exception:
            return False
    
    def _update_plot(self):
        """Update the convergence plot with new data"""
        try:
            # Read CSV file
            df = read_csv_file(self.csv_file)
            
            if df is None or df.shape[1] < 2:
                return
            
            # Extract x and y data
            x = df.iloc[:, 0].values
            y = df.iloc[:, 1].values
            
            # Filter out invalid values
            valid_mask = np.isfinite(x) & np.isfinite(y)
            x = x[valid_mask]
            y = y[valid_mask]
            
            if len(x) == 0 or len(y) == 0:
                return
            
            # Update series data (check if widget exists)
            if dpg.does_item_exist(self.series_tag):
                dpg.set_value(self.series_tag, [x.tolist(), y.tolist()])
            
            # Update axis limits
            self._update_axis_limits(x, y)
        
        except Exception as e:
            # Silently fail during updates
            pass
    
    def _update_axis_limits(self, x: np.ndarray, y: np.ndarray):
        """
        Update axis limits based on data
        
        Args:
            x: X-axis data
            y: Y-axis data
        """
        try:
            if len(x) == 0 or len(y) == 0:
                return
            
            # X-axis limits
            x_min = 0
            x_max = np.max(x) + 10
            
            if dpg.does_item_exist(self.x_axis_tag):
                dpg.set_axis_limits(self.x_axis_tag, x_min, x_max)
            
            # Y-axis limits (log scale)
            positive_y = y[y > 0]
            
            if len(positive_y) == 0:
                return
            
            y_min = np.nanmin(positive_y)
            y_max = np.nanmax(positive_y)
            
            # Add padding in log scale (one order of magnitude)
            y_min_log = np.floor(np.log10(y_min)) - 1
            y_max_log = np.ceil(np.log10(y_max)) + 0.2
            
            y_min = 10 ** y_min_log
            y_max = 10 ** y_max_log
            
            if dpg.does_item_exist(self.y_axis_tag):
                dpg.set_axis_limits(self.y_axis_tag, y_min, y_max)
        
        except Exception as e:
            if dpg.is_dearpygui_running():
                logger.debug(f"Error updating axis limits: {e}")
    
    def reset(self):
        """Reset the plot to initial state"""
        try:
            if dpg.does_item_exist(self.series_tag):
                dpg.set_value(self.series_tag, [[], []])
            
            if dpg.does_item_exist(self.x_axis_tag):
                dpg.set_axis_limits(self.x_axis_tag, 0, 100)
            
            if dpg.does_item_exist(self.y_axis_tag):
                dpg.set_axis_limits(self.y_axis_tag, 1e-10, 1)
            
            self.last_modified_time = 0
            
            logger.info("Convergence plot reset")
        
        except Exception as e:
            logger.log_exception(e, "Error resetting convergence plot")
    
    def force_update(self):
        """Force an immediate update of the plot"""
        try:
            self.last_modified_time = 0  # Force re-read
            self._update_plot()
        except Exception as e:
            logger.log_exception(e, "Error forcing convergence plot update")
    
    def is_running(self) -> bool:
        """
        Check if monitor is currently running
        
        Returns:
            True if monitoring, False otherwise
        """
        return self.is_monitoring
    
    def get_latest_data(self) -> Optional[tuple]:
        """
        Get the latest convergence data
        
        Returns:
            Tuple of (x_data, y_data) or None if no data available
        """
        try:
            if not os.path.exists(self.csv_file):
                return None
            
            df = read_csv_file(self.csv_file)
            
            if df is None or df.shape[1] < 2:
                return None
            
            x = df.iloc[:, 0].values
            y = df.iloc[:, 1].values
            
            return (x, y)
        
        except Exception:
            return None
    
    def get_convergence_status(self) -> dict:
        """
        Get current convergence status information
        
        Returns:
            Dictionary with convergence metrics
        """
        data = self.get_latest_data()
        
        if data is None:
            return {
                "has_data": False,
                "iterations": 0,
                "current_error": None,
                "min_error": None,
            }
        
        x, y = data
        
        if len(y) == 0:
            return {
                "has_data": False,
                "iterations": 0,
                "current_error": None,
                "min_error": None,
            }
        
        return {
            "has_data": True,
            "iterations": int(x[-1]) if len(x) > 0 else 0,
            "current_error": float(y[-1]) if len(y) > 0 else None,
            "min_error": float(np.min(y[y > 0])) if np.any(y > 0) else None,
            "is_converging": self._is_converging(y),
        }
    
    def _is_converging(self, y: np.ndarray, window: int = 10) -> bool:
        """
        Check if the solution is converging
        
        Args:
            y: Error values
            window: Window size for trend analysis
        
        Returns:
            True if converging, False otherwise
        """
        if len(y) < window:
            return False
        
        try:
            # Check if last 'window' values show decreasing trend
            recent = y[-window:]
            return recent[-1] < recent[0]
        except Exception:
            return False
    
    def cleanup(self):
        """Clean up resources"""
        self.stop()


# Global convergence monitor instance
convergence_monitor = ConvergenceMonitor()