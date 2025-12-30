"""
Font management utilities for MeMPhyS GUI

Handles font discovery, loading, and management across different platforms.
"""

import platform
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import dearpygui.dearpygui as dpg

from src.config import FONT_SIZES, FONT_PREFERENCES, FONT_SEARCH_PATHS
from src.core import app_state, logger


def find_system_font(preferred: Optional[List[str]] = None) -> str:
    """
    Find a suitable system font across different platforms
    
    Args:
        preferred: List of preferred font filenames to search for first
    
    Returns:
        Path to a valid font file
    
    Raises:
        FileNotFoundError: If no suitable font is found
    """
    system = platform.system()
    
    # Get search paths for this platform
    search_paths = FONT_SEARCH_PATHS.get(system, [])
    
    # Get default fonts for this platform
    defaults = FONT_PREFERENCES.get(system, [])
    
    # Combine preferred and defaults
    candidates = (preferred or []) + defaults
    
    logger.debug(f"Searching for fonts on {system}")
    logger.debug(f"Candidates: {candidates}")
    
    # Search in platform-specific directories
    for folder in search_paths:
        if not folder.exists():
            continue
            
        for font_name in candidates:
            font_path = folder / font_name
            if font_path.exists():
                logger.success(f"Found font: {font_path}")
                return str(font_path)
    
    # Fallback: search for any TTF file (Linux)
    if system == "Linux":
        logger.warning("No preferred fonts found, searching for any .ttf file")
        for folder in search_paths:
            if folder.exists():
                for font_file in folder.rglob("*.ttf"):
                    logger.info(f"Using fallback font: {font_file}")
                    return str(font_file)
    
    # If we get here, no font was found
    error_msg = f"No system font found. Tried: {candidates}"
    logger.error(error_msg)
    raise FileNotFoundError(error_msg)


def load_font(font_path: str, size: int) -> int:
    """
    Load a font with DearPyGUI
    
    Args:
        font_path: Path to the font file
        size: Font size in pixels
    
    Returns:
        Font ID from DearPyGUI
    
    Raises:
        FileNotFoundError: If font file doesn't exist
        RuntimeError: If font loading fails
    """
    if not Path(font_path).exists():
        raise FileNotFoundError(f"Font file not found: {font_path}")
    
    try:
        with dpg.font_registry():
            font_id = dpg.add_font(font_path, size)
        
        logger.debug(f"Loaded font: {Path(font_path).name} @ {size}px (ID: {font_id})")
        return font_id
    
    except Exception as e:
        error_msg = f"Failed to load font {font_path} @ {size}px: {e}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)


def initialize_fonts(preferred_fonts: Optional[List[str]] = None) -> Dict[Tuple[str, int], int]:
    """
    Initialize all fonts for the application
    
    Loads fonts for all sizes defined in config and registers them
    with the application state.
    
    Args:
        preferred_fonts: Optional list of preferred font names
    
    Returns:
        Dictionary mapping (font_name, size) to font IDs
    """
    logger.info("Initializing fonts...")
    
    system = platform.system()
    fonts_to_load = preferred_fonts or FONT_PREFERENCES.get(system, [])
    
    font_registry = {}
    
    for font_name in fonts_to_load:
        try:
            # Find the font file
            font_path = find_system_font([font_name])
            
            # Load all sizes
            for size in FONT_SIZES:
                try:
                    font_id = load_font(font_path, size)
                    font_registry[(font_name, size)] = font_id
                    
                    # Register with app state
                    app_state.register_font(font_name, size, font_id)
                    
                except Exception as e:
                    logger.warning(f"Could not load {font_name} @ {size}px: {e}")
            
            logger.success(f"Loaded font family: {font_name}")
            
        except FileNotFoundError:
            logger.warning(f"Font not found: {font_name}, trying next...")
            continue
    
    if not font_registry:
        logger.error("No fonts could be loaded!")
        raise RuntimeError("Failed to load any fonts")
    
    logger.success(f"Initialized {len(font_registry)} font variants")
    return font_registry


def set_default_font(font_name: str, size: int) -> bool:
    """
    Set the default font for the application
    
    Args:
        font_name: Name of the font
        size: Font size
    
    Returns:
        True if successful, False otherwise
    """
    font_id = app_state.get_font_id(font_name, size)
    
    if font_id is None:
        logger.error(f"Font not found: {font_name} @ {size}px")
        return False
    
    try:
        dpg.bind_font(font_id)
        app_state.current_font_name = font_name
        app_state.current_font_size = size
        logger.success(f"Set default font: {font_name} @ {size}px")
        return True
    
    except Exception as e:
        logger.error(f"Failed to set default font: {e}")
        return False


def change_font(font_name: str, size: int) -> bool:
    """
    Change the current font
    
    Args:
        font_name: Name of the font to change to
        size: Font size to change to
    
    Returns:
        True if successful, False otherwise
    """
    # Check if font exists in registry
    font_id = app_state.get_font_id(font_name, size)
    
    if font_id is None:
        logger.warning(f"Font not loaded: {font_name} @ {size}px, attempting to load...")
        
        # Try to load the font
        try:
            font_path = find_system_font([font_name])
            font_id = load_font(font_path, size)
            app_state.register_font(font_name, size, font_id)
        except Exception as e:
            logger.error(f"Could not load font: {e}")
            return False
    
    # Apply the font
    return set_default_font(font_name, size)


def get_available_fonts() -> List[str]:
    """
    Get list of available font names in the current system
    
    Returns:
        List of font names (without extensions)
    """
    system = platform.system()
    return FONT_PREFERENCES.get(system, [])


def get_available_sizes() -> List[int]:
    """
    Get list of available font sizes
    
    Returns:
        List of font sizes
    """
    return FONT_SIZES.copy()


def get_current_font() -> Tuple[str, int]:
    """
    Get the current font name and size
    
    Returns:
        Tuple of (font_name, font_size)
    """
    return (app_state.current_font_name, app_state.current_font_size)


def list_system_fonts(search_paths: Optional[List[Path]] = None) -> List[Path]:
    """
    List all font files available on the system
    
    Args:
        search_paths: Optional list of paths to search. If None, uses platform defaults
    
    Returns:
        List of Paths to font files
    """
    if search_paths is None:
        system = platform.system()
        search_paths = FONT_SEARCH_PATHS.get(system, [])
    
    font_files = []
    font_extensions = [".ttf", ".ttc", ".otf"]
    
    for folder in search_paths:
        if not folder.exists():
            continue
        
        for ext in font_extensions:
            font_files.extend(folder.rglob(f"*{ext}"))
    
    return font_files


def validate_font_file(font_path: str) -> bool:
    """
    Validate that a font file exists and is readable
    
    Args:
        font_path: Path to the font file
    
    Returns:
        True if valid, False otherwise
    """
    try:
        path = Path(font_path)
        return path.exists() and path.is_file() and path.suffix.lower() in [".ttf", ".ttc", ".otf"]
    except Exception:
        return False


def get_font_info() -> Dict[str, any]:
    """
    Get information about the current font configuration
    
    Returns:
        Dictionary with font information
    """
    current_name, current_size = get_current_font()
    
    # Count registered fonts
    registered_count = len(app_state._font_registry) if hasattr(app_state, '_font_registry') else 0
    
    return {
        "current_font": current_name,
        "current_size": current_size,
        "registered_fonts": registered_count,
        "available_fonts": get_available_fonts(),
        "available_sizes": get_available_sizes(),
        "platform": platform.system(),
    }


def create_font_with_fallback(primary_font: str, size: int, 
                              fallback_fonts: Optional[List[str]] = None) -> Optional[int]:
    """
    Create a font with automatic fallback to alternatives
    
    Args:
        primary_font: Primary font to try
        size: Font size
        fallback_fonts: List of fallback fonts to try if primary fails
    
    Returns:
        Font ID or None if all attempts fail
    """
    fonts_to_try = [primary_font]
    if fallback_fonts:
        fonts_to_try.extend(fallback_fonts)
    
    # Add system defaults as final fallback
    system = platform.system()
    fonts_to_try.extend(FONT_PREFERENCES.get(system, []))
    
    for font_name in fonts_to_try:
        try:
            font_path = find_system_font([font_name])
            font_id = load_font(font_path, size)
            app_state.register_font(font_name, size, font_id)
            logger.info(f"Successfully loaded font: {font_name}")
            return font_id
        except Exception as e:
            logger.debug(f"Could not load {font_name}: {e}")
            continue
    
    logger.error("All font loading attempts failed")
    return None