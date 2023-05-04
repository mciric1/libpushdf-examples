#!/bin/bash

# Quick setup script for building the fgb native code for libpushdf
cd libpushdf && python3 fgb/setup.py build_libfgb && python3 fgb/setup.py build
