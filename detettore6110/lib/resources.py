#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Resource management for detettore6110.
Handles access to reference genomes and IS target sequences.
"""

import os


def get_resource_path(relative_path):
    """
    Get the absolute path to a resource file.
    
    Resources are included with the package, but users can also provide
    their own paths. This function returns the package resource path
    if it exists, otherwise returns the relative path as-is (allowing
    user-provided paths).
    
    Parameters
    ----------
    relative_path : str
        Relative path to the resource file (e.g., 'resources/reference/MTBC0_v1.1.fasta')
    
    Returns
    -------
    str
        Absolute path to the resource file
    """
    # Try package resource path first
    package_dir = os.path.dirname(os.path.abspath(__file__))
    package_resource_path = os.path.join(package_dir, relative_path)
    
    if os.path.exists(package_resource_path):
        return package_resource_path
    
    # If not found in package, assume it's a user-provided path
    # Return as-is (could be absolute or relative to cwd)
    return relative_path


def get_default_reference():
    """Get the default reference genome path (MTBC0 v1.1)."""
    return get_resource_path('resources/reference/MTBC0_v1.1.fasta')


def get_default_is_target():
    """Get the default IS target path (IS6110)."""
    return get_resource_path('resources/is_targets/IS6110.fasta')


def get_default_annotation():
    """Get the default annotation path (MTBC0 v1.1 PGAP annotation)."""
    return get_resource_path('resources/reference/MTBC0v1.1_PGAP_annot.gff')
