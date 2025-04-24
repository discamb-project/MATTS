#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <vector>
#include <sstream>
#include <iomanip>

#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/MATTS_BankReader.h"

using namespace std;
using namespace discamb;

/*
ENTRY
  ID
    C306b
  COMMENT
    in UBDB2018: C306 splitted KJ
  NOI
    72
  ATOM DESCRIPTORS
    C1    CONNECTED_TO  O2,C3,C4          PLANARITY +  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING - IN_4_MEMBER_RING -
    O2    CONNECTED_TO  C1                PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING - IN_4_MEMBER_RING -
    C3    CONNECTED_TO  C1,X,X            PLANARITY +  PLANAR_RING_WITH_PLANAR_ATOMS *        IN_3_MEMBER_RING * IN_4_MEMBER_RING *
    C4    CONNECTED_TO  C1,X,X            PLANARITY +  PLANAR_RING_WITH_PLANAR_ATOMS *        IN_3_MEMBER_RING * IN_4_MEMBER_RING *
  LOCAL COORDINATE SYSTEM
    Z O2 X C3 R
  SYMMETRY
    mm2
  PARAMETER MODIFICATION DATE
    Fri Dec  3 10:28:39 2021
  MULTIPOLE MODEL PARAMETERS
    PVAL          4.052(89) KAPPA        0.9954(71) KPRIM         0.871(20)
    PLMS   1  0     0.068(37) PLMS   2  0     0.206(29) PLMS   2  2     0.212(20) PLMS   3  0     0.321(28)
    PLMS   3  2    -0.235(18) PLMS   4  0     0.025(17) PLMS   4  2     0.042(19) PLMS   4  4     0.024(11)
*/

void packType(
    const vector<string> &typeLines,
    vector<string> &packed)
{
    /*
    vector<string> words;
    packed.clear();
    string line = typeLines[2];
    line = string_utilities::trim(line) + " "; // ID
    int idx = 3;
    while (typeLines[idx].find("ATOM DESCRIPTORS") == string::npos)
        idx++;
    idx++;
    while (typeLines[idx].find("LOCAL COORDINATE SYSTEM") == string::npos)
    {    
        // 0        1            2                     3     4      5                        6    7             8     9            10
        // C1    CONNECTED_TO  O2, C3, C4          PLANARITY + PLANAR_RING_WITH_PLANAR_ATOMS - IN_3_MEMBER_RING - IN_4_MEMBER_RING -
        string_utilities::split(typeLines[idx++], words);
        for (int i = 0; i < 11; i += 2)
            line += words[i] + " ";
    }
    idx++;
    line += "LCS " + string_utilities::trim(string(typeLines[idx++]));
    idx++;
    line += " " + string_utilities::trim(string(typeLines[idx++]));
    
    
    
    string_utilities::split(typeLines[idx++], words);
    if (words[0] == "CHIRALITY")
    {
        line += " CH ";
        line += " " + string_utilities::trim(string(typeLines[idx++]));
    }

    packed.push_back(line);
    line.clear();
    
    while (typeLines[idx].find("MULTIPOLE MODEL PARAMETERS") == string::npos)
        idx++;
    idx++;
    string_utilities::split(typeLines[idx], words);
    for (int i = 1; i < 6; i += 2)
        line += words[i] + " ";

    for (idx = idx + 1; idx < typeLines.size(); idx++)
    {
        //PLMS   3  2    -0.235(18) PLMS   4  0     0.025(17) PLMS   4  2     0.042(19) PLMS   4  4     0.024(11)
        string_utilities::split(typeLines[idx], words);
        int n = words.size() / 4;
        for (int i = 0; i < n; i++)
            line += words[4 * i + 1] + " " + words[4 * i + 2] + " " + words[4 * i + 3] + " ";
    }
    packed.push_back(line);
    */
}

void unpackType(const string &line1, const string& line2, vector<string> &unpacked)
{
    vector<string> words;
    unpacked.clear();

    string_utilities::split(line1, words);

    unpacked.push_back("ENTRY");
    unpacked.push_back("  ID");
    unpacked.push_back("  " + words[0]);
    unpacked.push_back("  COMMENT");
    unpacked.push_back("    -");
    unpacked.push_back("  ATOM DESCRIPTORS");
    auto iterator_lcs = find(words.begin(), words.end(), "LCS");
    if (iterator_lcs == words.end())
        on_error::throwException("invalid format of packed bank", __FILE__, __LINE__);

    int lcs_word_idx = distance(words.begin(), iterator_lcs);
    int nAtoms = (lcs_word_idx - 1) / 6;
    if((lcs_word_idx - 1) % 6 != 0)
        on_error::throwException("invalid format of packed bank", __FILE__, __LINE__);

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        stringstream ss;
        ss << left;
        ss << "    " << setw(6) << words[atomIdx * 6 + 1]
            << "CONNECTED_TO  " << setw(18) << words[atomIdx * 6 + 2]
            << "PLANARITY " << setw(3) << words[atomIdx * 6 + 3]
            << "PLANAR_RING_WITH_PLANAR_ATOMS " << setw(9) << words[atomIdx * 6 + 4]
            << "IN_3_MEMBER_RING " << setw(2) << words[atomIdx * 6 + 5]
            << "IN_4_MEMBER_RING " << words[atomIdx * 6 + 6];
        string line;
        getline(ss, line);
        unpacked.push_back(line);
    }
    unpacked.push_back("LOCAL COORDINATE SYSTEM");
    unpacked.push_back("    ");
    for (int i = 0; i < 5; i++)
        unpacked.back() += words[lcs_word_idx + 1 + i] + (i==4? "": " ");
    unpacked.push_back("  SYMMETRY");
    unpacked.push_back("    " + words[lcs_word_idx + 6]);
    if (words.size() > lcs_word_idx + 7)
    {
        unpacked.push_back("  CHIRALITY");
        vector<string> chirality_def(words.begin()+ lcs_word_idx + 8, words.end());
        unpacked.push_back("    " + string_utilities::merge(chirality_def, ' '));
    }
    string_utilities::split(line2, words);
    unpacked.push_back("  MULTIPOLE MODEL PARAMETERS");
    {
        unpacked.push_back("");
        stringstream ss;
        ss  << "    PVAL" << setw(18) << words[0]
            << " KAPPA" << setw(18) << words[1]
            << " KPRIM" << setw(18) << words[2];
        getline(ss, unpacked.back());
    }
    int nPlm = (words.size() - 3) / 3;
    vector<string> plmStrings(nPlm);
    for (int i = 0; i < nPlm; i++)
    {
        stringstream ss;
        ss << "PLMS " << setw(3) << words[3 + 3 * i] << setw(3) << words[3 + 3 * i + 1]
            << setw(14) << words[3 + 3 * i + 2] << " ";
        getline(ss, plmStrings[i]);
    }
    int nPlmLines = nPlm / 3;
    if (nPlm % 3 != 0)
        nPlmLines++;

    
    for (int i = 0; i < nPlmLines; i++)
    {
        string line = "    ";
        int nPlmInLine = min(3,nPlm - 3 * i);
        for (int j = 0; j < nPlmInLine; j++)
            line += plmStrings[3 * i + j];
        unpacked.push_back(line);
    }

}

void unpack(const vector<string>& packed, vector<string>& unpacked)
{
    int idx = 0;
    unpacked.clear();
    while (packed[idx] != "ATOM TYPES")
        unpacked.push_back(packed[idx++]);
    unpacked.push_back("");
    idx++;
    int nTypes = (packed.size() - idx) / 2;
    for (int i = 0; i < nTypes; i++)
    {
        vector<string> unpackedLines;
        unpackType(packed[idx + 2 * i], packed[idx + 2 * i + 1], unpackedLines);
        unpacked.insert(unpacked.end(), unpackedLines.begin(), unpackedLines.end());
        unpacked.push_back("");
    }

}

void pack(
    const string& fileName,
    vector<string> &packed)
{
    packed.clear();
    ifstream in(fileName);
    string line;
    vector<string> typeLines, otherLines, words, packedType;
    bool readsType = false;
    
    while (in.good())
    {
        getline(in, line);
        if (!line.empty())
            if (line[0] == '#')
                continue;
        if (line.empty())
            continue;
        string_utilities::split(line, words);

        if (words.size() == 1)
            if (words[0] == "ENTRY")
            {
                readsType = true;
                if (!typeLines.empty())
                {
                    packType(typeLines, packedType);
                    packed.insert(packed.end(), packedType.begin(), packedType.end());
                }
                typeLines.clear();
            }

        if (!words.empty())
        {
            if (!readsType)
                otherLines.push_back(line);
            else
                typeLines.push_back(line);
        }
    }
    if(!typeLines.empty())
    {
        packType(typeLines, packedType);
        packed.insert(packed.end(), packedType.begin(), packedType.end());
    }

    in.close();
    ofstream out("out");
    for (auto &line: otherLines)
        out << "\"" << line << "\",\n";
    out << "\"ATOM TYPES\",\n";
    for (int i = 0; i < packed.size(); i++)
        out << "\"" << packed[i] << "\",\n";
    out.close();
}

void readPacked(const string& fName, vector<string>& lines)
{
    lines.clear();
    ifstream in(fName);
    string line;
    while (getline(in, line))
    {
        if (line.size() > 2)
            line = line.substr(1, line.size() - 3);
        lines.push_back(line);
    }
    in.close();
}

bool compare(
    const vector<AtomType> &t1,
    const vector<AtomType>& t2, 
    const vector<AtomTypeHC_Parameters> &p1,
    const vector<AtomTypeHC_Parameters> &p2,
    const BankSettings & settings1,
    const BankSettings & settings2)
{

    // types
    if (t1.size() != t2.size())
        return false;
    int nTypes = t1.size();

    for (int i = 0; i < nTypes; i++)
    {
        //cout << t1[i].id << endl;
        cout << i << endl;
        int nAtoms = t1[i].atoms.size();
        if (nAtoms != t2[i].atoms.size())
            return false;

        for (int j = 0; j < nAtoms; j++)
        {
            if (t1[i].atoms[j].label != t2[i].atoms[j].label)
                return false;
            if (t1[i].atoms[j].atomic_number_range != t2[i].atoms[j].atomic_number_range)
                return false;
            if (t1[i].atoms[j].anyAtomicNumber != t2[i].atoms[j].anyAtomicNumber)
                return false;
            if (t1[i].atoms[j].planar != t2[i].atoms[j].planar)
                return false;
            if (t1[i].atoms[j].neighborsAtomicNumbers != t2[i].atoms[j].neighborsAtomicNumbers)
                return false;
            if (t1[i].atoms[j].neighborsAtomicNumberRanges != t2[i].atoms[j].neighborsAtomicNumberRanges)
                return false;
            if (t1[i].atoms[j].nNeighbours != t2[i].atoms[j].nNeighbours)
                return false;
            if (t1[i].atoms[j].ringInfo.in3Ring != t2[i].atoms[j].ringInfo.in3Ring)
                return false;
            if (t1[i].atoms[j].ringInfo.in4Ring != t2[i].atoms[j].ringInfo.in4Ring)
                return false;
            if (t1[i].atoms[j].ringInfo.inAnyAdditionalRing != t2[i].atoms[j].ringInfo.inAnyAdditionalRing)
                return false;
            if (t1[i].atoms[j].ringInfo.inRing != t2[i].atoms[j].ringInfo.inRing)
                return false;
            if (t1[i].atoms[j].ringInfo.labeledContainingRings != t2[i].atoms[j].ringInfo.labeledContainingRings)
                return false;
            if (t1[i].atoms[j].ringInfo.labeledNonContainingRings != t2[i].atoms[j].ringInfo.labeledNonContainingRings)
                return false;
            if (t1[i].atoms[j].ringInfo.n3rings != t2[i].atoms[j].ringInfo.n3rings)
                return false;
            if (t1[i].atoms[j].ringInfo.n4rings != t2[i].atoms[j].ringInfo.n4rings)
                return false;
            if (t1[i].atoms[j].ringInfo.nonLabeledContainingRings != t2[i].atoms[j].ringInfo.nonLabeledContainingRings)
                return false;
            if (t1[i].atoms[j].ringInfo.nonLabeledNonContainingRings != t2[i].atoms[j].ringInfo.nonLabeledNonContainingRings)
                return false;
            if (t1[i].atoms[j].fixedNumberOfNeighbors != t2[i].atoms[j].fixedNumberOfNeighbors)
                return false;

        }


        if (t1[i].connectivity != t2[i].connectivity)
            return false;
        if (t1[i].ringLabels != t2[i].ringLabels)
            return false;
        if (t1[i].ringSizes != t2[i].ringSizes)
            return false;
        if (t1[i].symmetry != t2[i].symmetry)
            return false;
        if (t1[i].id != t2[i].id)
            return false;
        if (t1[i].ringLabels != t2[i].ringLabels)
            return false;
        if (t1[i].chirality != t2[i].chirality)
            return false;
        if (t1[i].localCoordinateSystem.isR != t2[i].localCoordinateSystem.isR)
            return false;
        if (t1[i].localCoordinateSystem.lcs_axis_1_definition != t2[i].localCoordinateSystem.lcs_axis_1_definition)
            return false;
        if (t1[i].localCoordinateSystem.lcs_axis_2_definition != t2[i].localCoordinateSystem.lcs_axis_2_definition)
            return false;
        if (t1[i].localCoordinateSystem.lcs_axis_type_1 != t2[i].localCoordinateSystem.lcs_axis_type_1)
            return false;
        if (t1[i].localCoordinateSystem.lcs_axis_type_2 != t2[i].localCoordinateSystem.lcs_axis_type_2)
            return false;
        if (t1[i].localCoordinateSystem.lcs_coordinate_1 != t2[i].localCoordinateSystem.lcs_coordinate_1)
            return false;
        if (t1[i].localCoordinateSystem.lcs_coordinate_2 != t2[i].localCoordinateSystem.lcs_coordinate_2)
            return false;
    }

    if (p1.size() != p1.size())
        return false;

    for (int i = 0; i < p1.size(); i++)
    {
        cout << i << " " << t1[i].id << endl;
        if (p1[i].kappa != p2[i].kappa)
            return false;
        if (p1[i].kappa_prime != p2[i].kappa_prime)
            return false;
        if (p1[i].kappa_prime_sigma != p2[i].kappa_prime_sigma)
            return false;
        if (p1[i].p_lms != p2[i].p_lms)
            return false;
        if (p1[i].p_lms_sigma != p2[i].p_lms_sigma)
            return false;
        if (p1[i].p_lm_indices != p2[i].p_lm_indices)
            return false;
        if (p1[i].p_val != p2[i].p_val)
            return false;
        if (p1[i].p_val_sigma != p2[i].p_val_sigma)
            return false;
        if (p1[i].symmetry != p2[i].symmetry)
            return false;
    }
    if (settings1.descriptorsSettings.atomInRingMaxNeighbourCount != settings2.descriptorsSettings.atomInRingMaxNeighbourCount)
        return false;
    if (settings1.descriptorsSettings.atomInRingPlanarityThreshold != settings2.descriptorsSettings.atomInRingPlanarityThreshold)
        return false;
    if (settings1.descriptorsSettings.atomPlanarityThreshold != settings2.descriptorsSettings.atomPlanarityThreshold)
        return false;
    if (settings1.descriptorsSettings.covalentBondThreshold != settings2.descriptorsSettings.covalentBondThreshold)
        return false;
    if (settings1.descriptorsSettings.maxCcDistanceAromaticRing != settings2.descriptorsSettings.maxCcDistanceAromaticRing)
        return false;
    if (settings1.descriptorsSettings.maxCnDistanceAromaticRing != settings2.descriptorsSettings.maxCnDistanceAromaticRing)
        return false;
    if (settings1.descriptorsSettings.maxPlanarRing != settings2.descriptorsSettings.maxPlanarRing)
        return false;
    if (settings1.descriptorsSettings.ringPlanarityThreshold != settings2.descriptorsSettings.ringPlanarityThreshold)
        return false;

    if (settings1.min_n_instaces != settings2.min_n_instaces)
        return false;
    if (settings1.min_plm != settings2.min_plm)
        return false;
    if (settings1.nSigma != settings2.nSigma)
        return false;
    return true;
}

int main(int argc, char *argv[])
{

    try 
    {
        /*
        if (argc != 2)
            cout<< "ERROR: expected bank file name as an argument, output will be placed in file 'output'\n";
        vector<string> packed, packedFromFile, unpacked;
        pack(argv[1], packed);

        readPacked("out", packedFromFile);
        unpack(packedFromFile, unpacked);

        ofstream out_unpacked("out2");
        for (string& line : unpacked)
            out_unpacked << line << "\n";
        out_unpacked.close();

        MATTS_BankReader reader;
        vector<AtomType> atomTypes1, atomTypes2;
        vector<AtomTypeHC_Parameters> param1, param2;
        BankSettings settings1, settings2;
        reader.read(argv[1], atomTypes1, param1, settings1);
        reader.read("out2", atomTypes2, param2, settings2);
        if(compare(atomTypes1, atomTypes2, param1, param2, settings1, settings2))
            cout<< "packed and unpacked bank the same as original\n";
        else
            cout << "packed and unpacked bank differ from the original\n";
        return 0;
        */
		ofstream out("out");
		ifstream in(argv[1]);
		string line;
		while (in.good())
		{
			getline(in, line);
			if (!line.empty())
				if (line[0] == '#')
					continue;
			out << "\"" << line << "\",\n";
		}
		in.close();
		out.close();

    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}
