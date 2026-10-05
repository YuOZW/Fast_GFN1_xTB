/* Regression for the highest-axis fallback used by Hessian thermochemistry.
 * Build with GNU fortification to catch the original two-byte buffer overflow:
 * gcc -O2 -D_FORTIFY_SOURCE=2 tests/portability/test_symmetry_fallback.c -lm -o test
 */
#include <assert.h>
#include "../../symmetry/symmetry_i.c"

int main(void)
{
    int counts[21] = {0};
    NormalAxesCounts = counts;
    MaxAxisOrder = 20;
    NormalAxesCount = 2;
    counts[2] = counts[4] = 1;
    MaxRotAxis[0] = '\0';
    report_symmetry_elements_brief_Conly();
    assert(strcmp(MaxRotAxis, "C4") == 0);

    counts[20] = 1;
    MaxRotAxis[0] = '\0';
    report_symmetry_elements_brief_Conly();
    assert(strcmp(MaxRotAxis, "C20") == 0);

    memset(counts, 0, sizeof(counts));
    NormalAxesCount = 0;
    MaxRotAxis[0] = '\0';
    report_symmetry_elements_brief_Conly();
    assert(MaxRotAxis[0] == '\0');
    return 0;
}
