#include "MATTS_Default.h"

#include <vector>

using namespace std;



namespace discamb{
	void default_ubdb_bank_string(std::string &s)
	{
		vector<string> v = {
"",
"SETTINGS",
"  covalent bond threshold 0.4",
"  atom planarity threshold 0.1",
"  ring planarity threshold 0.1",
"  atom in ring planarity threshold 0.1",
"  atom in planar ring max number of neighbours 3",
"  min plm 0.002",
"  n sigma 1.0",
"  minimal number of type instances 2",
"",
"",
"ENTRY",
"  ID",
"    H101",
"  COMMENT",
"    in UBDB2018: H101 KJ  ",
"  NOI",
"    10089",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2                PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,H,H,!H         PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 X !H(C2) R",
"  SYMMETRY",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Dec  3 10:28:39 2021",
"  MULTIPOLE MODEL PARAMETERS",
"    PVAL          1.092(38) KAPPA         1.093(13) KPRIM         1.149(19)",
"    PLMS   1  0     0.186(10) PLMS   2  0    0.0857(87)",
"",
"ENTRY",
"  ID",
"    H102",
"  COMMENT",
"    in UBDB2018: H102 KJ  ",
"  NOI",
"    7124",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2                PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,H,!H,!H        PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS *        IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 X !H(C2) R",
"  SYMMETRY",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Dec  3 10:28:39 2021",
"  MULTIPOLE MODEL PARAMETERS",
"    PVAL          1.064(45) KAPPA         1.098(16) KPRIM         1.140(31)",
"    PLMS   1  0     0.180(15) PLMS   2  0     0.086(11)",
"",
""};

	s.clear();
	for (auto const& line : v)
		s += line + string("\n");
	}
}
