"""
The top-level augur command which dispatches to subcommands.
"""
command_strings = [
	"gRNA_mutation",
	"editing_effect"]

import argparse
import girafr.gRNA_mutation
import girafr.editing_effect
import sys


parser = argparse.ArgumentParser(description='GiRAFR.')
parser.add_argument('command', default = "gRNA_mutation", choices = ["gRNA_mutation","editing_effect"], help = "command gRNA_mutation or editing_effect (default: %(default)s)")
parser.add_argument('--config', '-c', default="ConfigFile", help="configuration file (default: %(default)s)")
args = parser.parse_args()

def run(args):
	if args.command == 'gRNA_mutation':
		girafr.gRNA_mutation.run(args.config)
	else:
		girafr.editing_effect.run(args.config)

run(args)
