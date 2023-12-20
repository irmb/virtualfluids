r"""
=======================================================================================
 ____          ____    __    ______     __________   __      __       __        __
 \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
      \    \  |    |   ________________________________________________________________
       \    \ |    |  |  ______________________________________________________________|
        \    \|    |  |  |         __          __     __     __     ______      _______
         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/

  This file is part of VirtualFluids. VirtualFluids is free software: you can
  redistribute it and/or modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  SPDX-License-Identifier: GPL-3.0-or-later
  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

! \author Henry Korb
=======================================================================================
"""
from __future__ import annotations

from typing import overload

class ConfigurationFile:
    def __init__(self) -> None: ...
    def contains(self, key: str) -> bool: ...
    @overload
    def get_bool_value(self, key: str) -> bool: ...
    @overload
    def get_bool_value(self, key: str, default_value: bool) -> bool: ...
    @overload
    def get_double_value(self, key: str) -> float: ...
    @overload
    def get_double_value(self, key: str, default_value: float) -> float: ...
    @overload
    def get_float_value(self, key: str) -> float: ...
    @overload
    def get_float_value(self, key: str, default_value: float) -> float: ...
    @overload
    def get_int_value(self, key: str) -> int: ...
    @overload
    def get_int_value(self, key: str, default_value: int) -> int: ...
    @overload
    def get_string_value(self, key: str) -> str: ...
    @overload
    def get_string_value(self, key: str, default_value: str) -> str: ...
    @overload
    def get_uint_value(self, key: str) -> int: ...
    @overload
    def get_uint_value(self, key: str, default_value: int) -> int: ...
    def load(self, file: str) -> bool: ...
