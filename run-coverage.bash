#!/usr/bin/env bash

pytest --cov=semirings tests/test_semirings.py
coverage html
xdg-open htmlcov/index.html
