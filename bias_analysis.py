from parser_test import *

analysis = bias_finder()

analysis.quality_control()
analysis.link_lib()
analysis.print_link_lib()
