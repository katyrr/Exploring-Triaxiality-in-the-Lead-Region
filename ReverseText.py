#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 21:09:35 2024

@author: katyrr
"""

# read some text and overwrite with the line order backwards

filepath = "energy_levels.txt"


read_file = open(filepath, 'r')

lines = []

for l in read_file:
    lines.append(l.strip())
    
read_file.close()

lines.reverse()


backwards_text = ""

for l in lines:
    backwards_text += l
    backwards_text += "\n"
    print("hi")

print(backwards_text)

write_file = open(filepath, 'w')

write_file.write(backwards_text)