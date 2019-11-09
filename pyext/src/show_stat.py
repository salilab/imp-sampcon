from __future__ import print_function
import sys
import IMP.pmi.output


p = IMP.pmi.output.ProcessOutput(sys.argv[1])
fields = p.get_keys()
print("\n".join(sorted(fields)))
