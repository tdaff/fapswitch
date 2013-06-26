"""
Raw ansi logo that also inserts some version information.

"""

from . import __version__

LOGO = r""" [0;30;47;107m                          [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m                          [0m  [0;30;47;107m      __                                              [0m 
 [0;30;47;107m   [0;37;5;47;107m      [0;1;37;97;47m [0;37;5;40;100m [0;36;5;40;100m [0;1;30;90;47m    [0;1;37;97;47m [0;37;5;47;107m      [0;30;47;107m   [0m  [0;30;47;107m     / _| __ _ _ __  ___                              [0m
 [0;30;47;107m   [0;37;5;47;107m   [0;1;30;90;47m [0;30;5;40;100m [0;1;30;90;40m [0;34;40m [0;1;30;90;40m [0;30;5;40;100m   [0;1;30;90;40m [0;31;5;40;100m [0;30;5;40;100m [0;33;5;40;100m  [0;37;5;47;107m    [0;30;47;107m   [0m  [0;30;47;107m    | |_ / _` | '_ \/ __|     Fully automated         [0m 
 [0;30;47;107m   [0;37;5;47;107m  [0;36;5;40;100m [0;1;30;90;40m     [0;31;40m [0;30;5;40;100m [0;31;5;40;100m [0;33;5;40;100m [0;35;5;40;100m [0;31;5;40;100m  [0;30;5;40;100m [0;33;5;40;100m [0;37;5;47;107m   [0;30;47;107m   [0m  [0;30;47;107m    |  _| (_| | |_) \__ \     adsorption analysis     [0m 
 [0;30;47;107m   [0;37;5;47;107m [0;37;5;40;100m [0;30;5;40;100m [0;1;30;90;40m    [0;31;40m [0;33;5;40;100m [0;37;5;40;100m [0;1;37;97;47m     [0;1;33;93;47m [0;35;5;40;100m [0;37;5;47;107m   [0;30;47;107m   [0m  [0;30;47;107m    |_|  \__,_| .__/|___/     in porous solids        [0m 
 [0;30;47;107m   [0;37;5;47;107m [0;36;5;40;100m [0;31;40m [0;30;5;40;100m [0;31;5;40;100m   [0;33;5;40;100m [0;1;30;90;47m [0;37;5;47;107m  [0;1;33;93;47m  [0;1;37;97;47m  [0;37;5;47;107m [0;33;47m [0;37;5;47;107m   [0;30;47;107m   [0m  [0;30;47;107m              |_|                                     [0m 
 [0;30;47;107m   [0;37;5;47;107m [0;1;37;97;47m [0;1;30;90;40m [0;31;5;40;100m [0;1;37;97;47m  [0;1;30;90;47m  [0;37;5;47;107m      [0;1;33;93;5;47;107m [0;1;35;95;47m [0;37;5;47;107m    [0;30;47;107m   [0m  [0;30;47;107m    ================>>>>>     Version:                [0m 
 [0;30;47;107m   [0;37;5;47;107m  [0;37;5;40;100m [0;31;40m [0;33;5;40;100m [0;37;5;47;107m [0;1;37;97;47m [0;37;5;47;107m  [0;1;33;93;5;47;107m [0;37;5;47;107m    [0;1;35;95;5;41;101m [0;1;33;93;5;47;107m [0;37;5;47;107m    [0;30;47;107m   [0m  [0;30;47;107m    <<<<<================     XXX.XXX.XXX.XXXXX       [0m 
 [0;30;47;107m   [0;37;5;47;107m   [0;1;37;97;47m [0;1;33;93;47m [0;1;37;97;47m [0;1;33;93;5;47;107m [0;1;33;93;47m [0;37;5;45;105m [0;37;5;43;103m [0;1;35;95;5;47;107m [0;1;37;97;5;43;103m [0;1;35;95;5;47;107m [0;1;37;97;5;43;103m [0;37;5;47;107m [0;1;37;97;47m [0;37;5;47;107m    [0;30;47;107m   [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m   [0;37;5;47;107m [0;1;37;97;47m [0;35;47m [0;37;5;41;101m [0;1;31;91;47m  [0;1;37;97;5;45;105m [0;1;37;97;5;43;103m [0;37;5;47;107m [0;37;5;41;101m [0;1;31;91;47m [0;37;5;41;101m  [0;1;30;90;47m [0;1;37;97;47m [0;37;5;47;107m     [0;30;47;107m   [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m   [0;1;30;90;47m [0;35;5;41;101m [0;37;5;41;101m [0;1;35;95;5;41;101m [0;37;5;41;101m [0;1;35;95;5;41;101m [0;1;33;93;5;41;101m [0;1;35;95;5;41;101m  [0;1;33;93;47m [0;1;37;97;47m [0;1;33;93;47m [0;1;37;97;47m [0;37;5;47;107m       [0;30;47;107m   [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m   [0;35;5;41;101m [0;35;47m [0;35;5;41;101m [0;1;31;91;47m [0;1;35;95;5;41;101m [0;37;5;41;101m [0;1;35;95;5;41;101m [0;37;5;41;101m  [0;1;35;95;5;41;101m [0;1;33;93;5;41;101m [0;37;5;47;107m         [0;30;47;107m   [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m                          [0m  [0;30;47;107m                                                      [0m 
 [0;30;47;107m (c) Tom Daff             [0m  [0;30;47;107m                                                      [0m """

formatted_version = "{0:<17}".format(__version__)

LOGO = LOGO.replace('XXX.XXX.XXX.XXXXX', formatted_version)

