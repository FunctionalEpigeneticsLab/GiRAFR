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
	if len(argv) != 1:
		print('Usage: girafr gRNA_mutation or girafr editing_effect')
	if argv[0] not in command_strings:
		print('Unknown command', argv[0])
		print('Usage: girafr gRNA_mutation or girafr editing_effect')
		sys.exit(2)
	if argv[0] == 'gRNA_mutation':
		girafr.gRNA_mutation.run()
	else:
		girafr.editing_effect.run()

