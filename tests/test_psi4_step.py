#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psi4_step` package."""

import pytest
import psi4_step
import re  # noqa: F401


@pytest.fixture
def instance():
    instance = psi4_step.Psi4()
    instance._id = (1,)
    return instance


def test_construction():
    """Simplest test that we can make a Psi4 object"""
    instance = psi4_step.Psi4()
    assert (str(type(instance)) == "<class 'psi4_step.psi4.Psi4'>")


def test_version():
    """Test that the object returns a version"""
    instance = psi4_step.Psi4()
    result = instance.version
    assert isinstance(result, str) and len(result) > 0


def test_git_revision():
    """Test that the object returns a git revision"""
    instance = psi4_step.Psi4()
    result = instance.git_revision
    assert isinstance(result, str) and len(result) > 0


# def test_description_text_default(instance):
#     """Test the default description text"""

#     print(instance.description_text())
#     assert re.fullmatch(
#         (
#             r'Step 1: Psi4  [-+.0-9a-z]+\n'
#             r'    Create a 2 x 2 x 2 supercell from the current cell'
#         ), instance.description_text()
#     ) is not None
