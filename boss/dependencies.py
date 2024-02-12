import sys
import subprocess
from shutil import which


class Dependencies:

    def __init__(self):
        """
        Initialize the Dependencies class.

        :raises SystemExit: If any of the dependencies are not found in the system path.
        """
        self.dependencies = ["minimap2", "paf2gfa", "miniasm"]
        self.check_dependencies()


    def check_dependencies(self) -> None:
        """
        Check if all dependencies are present in the system path.

        :raises SystemExit: If any of the dependencies are not found in the system path.
        """
        for dependency in self.dependencies:
            setattr(self, dependency, find_exe(dependency))
            # logging.info(f'{dependency}: {getattr(self, dependency)}')
            if not getattr(self, dependency):
                sys.exit(f"Dependency {dependency} not found in path")



def find_exe(name: str) -> str | None:
    """
    Find the executable file with the specified name.

    :param name: The name of the executable.
    :return: The path to the executable, or None if not found.
    """
    # shutil.which seems to work mostly but is still not completely portable
    exe = which(name, path='/'.join(sys.executable.split('/')[:-1]))
    if not exe:
        exe = which(name)
    if not exe:
        exe = subprocess.run(f'which {name}', shell=True, capture_output=True, universal_newlines=True).stdout
    if not exe:
        return
    return exe.strip()

