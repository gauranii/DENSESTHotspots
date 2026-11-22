# -----------------------------------------------------------------------------
# File:    python/tests/test_every_function.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Every module of python/zinull has its test file, and every public
#          function of a module is called in that file.
# Usage:   python -m pytest python/tests/test_every_function.py
# -----------------------------------------------------------------------------
import ast
import os
import re

import pytest

from conftest import TESTS

PACKAGE = os.path.join(TESTS, "..", "zinull")
MODULES = sorted(file[:-3] for file in os.listdir(PACKAGE)
                 if file.endswith(".py") and file != "__init__.py")


def public_functions(module):
    """Names of the top-level functions of a module that do not start with
    an underscore (also those defined under `if`, as rmath.muladd)."""
    with open(os.path.join(PACKAGE, module + ".py")) as f:
        tree = ast.parse(f.read())
    nodes = list(tree.body)
    for node in tree.body:
        if isinstance(node, ast.If):
            nodes += node.body + node.orelse
    return sorted({node.name for node in nodes
                   if isinstance(node, ast.FunctionDef) and not node.name.startswith("_")})


def test_the_package_has_the_modules_of_r():
    r_modules = sorted(file[:-2] for file in os.listdir(os.path.join(TESTS, "..", "..", "R")))
    assert set(r_modules) <= set(MODULES)
    assert set(MODULES) - set(r_modules) == {"rng", "rmath", "rbase"}


@pytest.mark.parametrize("module", MODULES)
def test_every_public_function_is_tested(module):
    with open(os.path.join(TESTS, f"test_{module}.py")) as f:
        tests = f.read()
    functions = public_functions(module)
    assert functions
    missing = [name for name in functions if not re.search(rf"\b{name}\(", tests)]
    assert not missing, f"not called in test_{module}.py: {missing}"
