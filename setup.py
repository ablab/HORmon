# Copyright (c) 2021 Saint Petersburg State University
# This code is part of HORmon software
#
# HORmon is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 2,
# dated June 1991, as published by the Free Software Foundation.
#
# HORmon is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="HORmon",
    version="1.0",
    description="Centromere annotation tool",
    packages=find_packages(include=['HORmon', 'HORmon.*']),
    entry_points={
        'console_scripts': ['HORmon=HORmon.HORmon:main', 'monomer_inference=HORmon.monomer_inference:main']
    }
)
