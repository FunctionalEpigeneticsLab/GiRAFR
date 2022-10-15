"""
Stub function and module used as a setuptools entry point.
"""
import girafr
from sys import argv, exit

# Entry point for setuptools-installed script and bin/girafr dev wrapper.
def main():
	return girafr.run(argv[1:])

# Run when called as `python -m girafr`, here for good measure.
if __name__ == "__main__":
	exit( main() )
