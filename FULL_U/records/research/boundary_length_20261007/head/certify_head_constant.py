"""Use the existing negative-constant threshold cap with centered interiors.

This adds only a choice of already reviewed methods. The frozen original
driver and every shared evaluator remain unchanged. A constant-margin cap
does not require a negative quadratic Taylor coefficient.
"""
from pathlib import Path
import hashlib
import certify_head as driver


class FiniteConstantCentered(driver.FiniteCentered):
    def threshold_cap(self,b):
        return driver.FiniteSpin.threshold_cap(self,b)


driver.MODELS['constant_centered']=FiniteConstantCentered
original_hashes=driver.code_hashes


def code_hashes():
    result=original_hashes()
    result[str(Path(__file__).relative_to(driver.ROOT))]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result


driver.code_hashes=code_hashes

if __name__=='__main__':driver.main()
