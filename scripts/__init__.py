"""
The top-level augur command which dispatches to subcommands.
"""
command_strings = [
	"gRNA_mutation",
	"editing_effect"]

#COMMANDS = [importlib.import_module('girafr.process_gRNA_mutation', 'girafr.process_detect_editing_effect')]
import girafr.gRNA_mutation
import girafr.editing_effect
import sys

def run(argv):
	if argv not in command_strings:
		print('Unknown command', argv)
		sys.exit(2)
	if argv == 'gRNA_mutation':
		girafr.gRNA_mutation.run()
	else:
		girafr.editing_effect.run()

