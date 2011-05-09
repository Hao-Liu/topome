/*
 * topome - a molecular dynamics and topology analysis package for 
 * coordination systems
 *
 * /__ __X  _ \/  __\/  _ \/ \__/|/  __/ Copyright (C) 2010 - Liu Hao
 *   / \ | / \||  \/|| / \|| |\/|||  \  
 *   | | | \_/||  __/| \_/|| |  |||  /_  topome.c
 *   \_/ \____/\_/   \____/\_/  \|\____\
 *
 * This file is the main entrance of topome.
 * 
 *
 * topome is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * topome is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with topome; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, 
 * Boston, MA  02110-1301  USA
 */

#include <stdio.h>
#include <string.h>
#include <libintl.h>
#include <getopt.h>
#include "verbose.h"
#include "graphic.h"
#include "parse.h"
#include "forcefield.h"


int 
main (int argc, char **argv)
{
  // default to graphic mode
  int graphic_enabled=0, verbose_enabled=0;
    
  // system structure for global variables
  System system;

  // OpenGL Window handle
  GLWindow gl_window;
  
  char *input_file;
  
  // parse input options
  parse_option (&input_file, &graphic_enabled, &verbose_enabled,
                      argc, argv);
  
  memset (&system.forcefield, 0, sizeof(ForceField));
  init_forcefield ("frc/cvff.frc", &system.forcefield);
  // initialize system
  init_system (&system, input_file);

  // relax to min energy mode
	relax_system (&system);

  // init graphic Xwindow if graphic mode is on
  if (graphic_enabled)
  {
    init_graphic (&system, &gl_window);
  }

  // iterate simulation
  run_system (&system, verbose_enabled, graphic_enabled, &gl_window);

  // clean up 
  release_system (&system);
  if (graphic_enabled)
  {
    release_graphic (&gl_window);
  }
  
  return 0;
}
