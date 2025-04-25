#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/BasicUtilities/discamb_version.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/IO/xd_io.h"
#include "discamb/IO/shelx_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/structure_io.h"
//#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
#include "MATTS_Default.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/BasicUtilities/discamb_env.h"

#include "json.hpp"

#include <iostream>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <set>
#include <filesystem>

using namespace discamb;
using namespace std;


void add_to_discamb_log(const string& fileName, const string& message)
{
    ofstream out(fileName, ostream::app);
    out << message;
    out.close();
}


void makeWholeHklSet(
    const vector<Vector3i>& hkl0, 
    const vector<Matrix3i> &rotations, 
    const Matrix3i &twinMatrix,
    int nTwinComponents,
    vector<Vector3i> &hkl)
{
	hkl.clear();
	set<Vector3i> uniqueHkl;
    vector<Vector3i> hkl1 = hkl0;

    for (auto const& h : hkl0)
    {
        Vector3i _h = h;
        for (int i = 1; i < nTwinComponents; i++)
        {
            _h = twinMatrix * _h;
            hkl1.push_back(_h);
        }
    }

	for (auto const& rotation : rotations)
        for (auto const& h : hkl1)
            uniqueHkl.insert(h * rotation);

	hkl.assign(uniqueHkl.begin(), uniqueHkl.end());

}

void json_update(nlohmann::json& a, nlohmann::json& b)
{
    cout << "merging JSON files\n";

    for (auto it = b.begin(); it != b.end(); it++)
    {
        if (a.find(it.key()) == a.end())
            a[it.key()] = it.value();
        else
        {
            if (it->is_object())
                json_update(*a.find(it.key()), *it);
        }
    }
}

void defaultInput(
    bool electronScatteringIfNoAspherJsonFile,
    bool unitCellChargeIfNoAspherJsonFile,
    double unitCellCharge,
    nlohmann::json& jsonData,
    const std::string & jobName)
{
    nlohmann::json data;

    stringstream ss;
    ss <<
        "{\n"
        "    \"assignment info\" : \"print_to_discamb2tsc_log_file\"\n,"
        "    \"multipole cif\" : \"" + jobName + ".cif_rho\"";
    if (electronScatteringIfNoAspherJsonFile)
        ss << ",\n   \"electron_scattering\": true";

    if (unitCellChargeIfNoAspherJsonFile)
        ss << ",\n   \"unit cell charge\": " << unitCellCharge;

    ss << "\n"
          "}\n";

    ss >> data;

    jsonData = data;
}


std::shared_ptr<HcAtomBankStructureFactorCalculator> sfCalculatorFromJsonFile(
    const Crystal &crystal,
    bool &aspherJsonPresent,
    bool electronScatteringIfNoAspherJsonFile,
    bool unitCellChargeIfNoAspherJsonFile,
    double unitCellCharge,
    bool &electronScattering,
    const string &jobName)
{
    
    aspherJsonPresent = filesystem::exists(filesystem::path("aspher.json"));
    nlohmann::json jsonData;
    
    if (aspherJsonPresent)
    {
        ifstream jsonFileStream("aspher.json");
        
        if (jsonFileStream.good())
        {
            nlohmann::json jsonDataFromFile;
            jsonFileStream >> jsonDataFromFile;
            auto data = jsonDataFromFile.find("form factor engine");
            if (data != jsonDataFromFile.end())
                jsonData = jsonDataFromFile["form factor engine"]["data"];
            else
                jsonData = jsonDataFromFile;
        }
        else
        {
            on_error::throwException("can not read aspher.json file, expected to be present in the current directory", __FILE__, __LINE__);
            return std::shared_ptr<HcAtomBankStructureFactorCalculator>(nullptr);
        }
    }
    else
        defaultInput(electronScatteringIfNoAspherJsonFile, unitCellChargeIfNoAspherJsonFile, unitCellCharge, jsonData, jobName);


    
    // electrons? bank path give?
    electronScattering = false;
    bool hasBankPath = false;

    electronScattering = jsonData.value("electron_scattering", electronScattering);
    electronScattering = jsonData.value("electron scattering", electronScattering);
    
    string radiation = electronScattering ? "electron" : "X-ray";
    clog << "Atomic form factor for " << radiation << " scattering calculated with\nHansen - Coppens model parameterized with MATTS databank.\n\n";

    if (jsonData.find("bank path") != jsonData.end())
        hasBankPath = true;

    if(hasBankPath)
        return make_shared<HcAtomBankStructureFactorCalculator>(crystal, jsonData);

    string bankString;
    default_ubdb_bank_string(bankString);
    return make_shared<HcAtomBankStructureFactorCalculator>(crystal, jsonData, bankString);
}

    bool hklShouldBeTakenFromTsc(
        string& tscFile)
    {
        nlohmann::json data;

        if (!filesystem::exists(filesystem::path("aspher.json")))
            return false;

        ifstream jsonFileStream("aspher.json");
        if (!jsonFileStream.good())
            on_error::throwException("cannot read aspher.json", __FILE__, __LINE__);

        if (jsonFileStream.good())
            jsonFileStream >> data;

        jsonFileStream.close();
        tscFile = data.value("hkl from tsc", string());
        if (tscFile.empty())
            return false;
        return true;
    }



void getStructureAndHklFile(
    const vector<string> &arguments,
    string &structureFile,
    string &hklFile,
    bool &cifAndHklFromDirectorySearch)
{
    // structure file suffix - cif, res ins dis
    // hkl suffix - hkl, fcf, tsc

    string tscFile;
    bool hklFromTsc = hklShouldBeTakenFromTsc(tscFile);

    if (arguments.size() == 2)
    {
        structureFile = arguments[0];
        hklFile = arguments[1];
        if (hklFromTsc)
            hklFile = tscFile;
        cifAndHklFromDirectorySearch = false;
        return;
    }

    if (hklFromTsc && arguments.size() == 1)
    {
        structureFile = arguments[0];
        hklFile = tscFile;
        cifAndHklFromDirectorySearch = false;
        return;
    }

    cifAndHklFromDirectorySearch = true;

    structureFile = file_system_utilities::find_newest_file("cif");

    if(hklFromTsc)
        hklFile = tscFile;
    else
        hklFile = file_system_utilities::find_newest_file("hkl");


    if (structureFile.empty())
        on_error::throwException("CIF file not found", __FILE__, __LINE__);
    if (hklFile.empty())
        on_error::throwException("hkl file not found", __FILE__, __LINE__);
    
}


void addMattsCitation(string mattsVersion)
{

    clog<< "Citation\n"
        << "If you used discambMATTS2tsc in your work please add the following text\n"
        << "to the resulting.cif file :\n"
        << "_refine_special_details\n"
        << ";\n"
        << "TAAM / MATTS refinement.\n"
        << "Uses aspherical atomic scattering factors computed by the DiSCaMB library\n"
        << "(Chodkiewicz et.al., J.Appl.Cryst., 2018, 51, 193 - 199)\n"
        << "from multipole model\n"
        << "(Hansen & Coppens, Acta Cryst.A, 1978, 34, 909 - 921)\n"
        << "parametrized using the MATTS2021 data bank\n"
        << "(Jha, et al., J.Chem.Inf.Model., 2022, 62, 3752 - 3765,\n"
        << "    Rybicka, et al., J.Chem.Inf.Model., 2022, 62, 3766 - 3783)\n"
        << "Refinement performed with.tsc generated by discamb2TAAMtsc " << mattsVersion << "\n"
        << ";\n"
        << "and in your publication, please, cite:\n"
        << "Chodkiewicz et.al., J.Appl.Cryst., 2018, 51, 193 - 199\n"
        << "Hansen & Coppens, Acta Cryst.A, 1978, 34, 909 - 921\n"
        << "Jha, et al., J.Chem.Inf.Model., 2022, 62, 3752 - 3765\n"
        << "Rybicka, et al., J.Chem.Inf.Model., 2022, 62, 3766 - 3783\n";
}

void print_bank()
{
    ofstream out("matts_data.txt");
    string bankContent;
    default_ubdb_bank_string(bankContent);
    out << bankContent << endl;
    out.close();
}

void read_structure(
    const string structureFile,
    Crystal& crystal,
    int& nTwinComponents,
    Matrix3i& twinMatrix)
{
    nTwinComponents = 1;
    string suffix = filesystem::path(structureFile).extension().string();
    if (suffix == string(".res") || suffix == string(".ins"))
    {
        map<string, string> data;
        shelx_io::read(structureFile, crystal, data);
        if (data.find("TWIN") != data.end())
        {
            vector<string> words;
            string_utilities::split(data.find("TWIN")->second, words);
            int nWords = words.size();
            if (nWords != 9 && nWords != 10)
                on_error::throwException("error when processing TWIN instruction in file '" + structureFile + "'", __FILE__, __LINE__);
            nWords == 9 ? nTwinComponents = 2 : nTwinComponents = stoi(words[9]);
            twinMatrix.set(stoi(words[0]), stoi(words[1]), stoi(words[2]),
                           stoi(words[3]), stoi(words[4]), stoi(words[5]),
                           stoi(words[6]), stoi(words[7]), stoi(words[8]));
        }
    }
    else
        structure_io::read_structure(structureFile, crystal);

}

void readHkl(
    const string& hklFile,
    vector<Vector3i>& hkls,
    bool shelxFreeFormat)
{
    hkls.clear();
    vector<int> batchNumbers;
    vector<double> intensities, sigmas;

    if (filesystem::path(hklFile).extension().string() == string(".tsc"))
    {
        vector<string> atomLabels;
        vector<vector<complex<double> > > ff;
        tsc_io::read_tsc(hklFile, atomLabels, hkls, ff);
        hkls.push_back(Vector3i(0, 0, 0));
    }
    else
        hkl_io::readShelxHkl(hklFile, hkls, intensities, sigmas, batchNumbers, shelxFreeFormat);

}

int main(int argc, char *argv[])
{
    Crystal crystal;
    UnitCellContent ucContent;
    UnitCellContent::AtomID atomID;
    Vector3d positionCartesian,positionFractional;
    vector<string> symbols;
    vector<vector<UnitCellContent::AtomID> > molecules;
    vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
    ofstream out;
    int nAtoms;
    vector<string> arguments; 
    vector<string> options;
    map<string, string> optionsWithValues;
    bool defaultTaamRun = false;
    bool electronScattering = false;


    try {
        try {

            ofstream logFile("discambMATTS2tsc.log");
            auto clog_orginal_rdbuf = std::clog.rdbuf();
            std::clog.rdbuf(logFile.rdbuf());


            string versionString = discamb_version::version();

            time_t time_now = std::chrono::system_clock::to_time_t(chrono::system_clock::now());
            string date = __DATE__;
            string time = __TIME__;
            string discambMATTS2tscVersion = "3.009";
            string header = string("\n  discambMATTS2tsc, version ") + discambMATTS2tscVersion + string("\n  compiled on ") +
                date + string(" , ") + time + string(".\n\n");

            cout << header;
            
            bool useSymmetrySection = false;
            bool newImplementation = true;


            // read crystal data

            string structureFile, hklFile;

            parse_cmd::get_args_and_options(argc, argv, arguments, options, optionsWithValues);
            
            if (find(options.begin(), options.end(), "-pb") != options.end())
            {
                print_bank();
                return 0;
            }
            
            bool cifAndHklFromDirectorySearch;
            getStructureAndHklFile(arguments, structureFile, hklFile, cifAndHklFromDirectorySearch);

            string jobName = filesystem::path(structureFile).stem().string();
            
            defaultTaamRun = arguments.empty() && cifAndHklFromDirectorySearch && (!filesystem::exists(filesystem::path("aspher.json")));

            filesystem::path structureFilePath(structureFile);
            bool shelxFreeFormat = false;
            bool electronScatteringIfNoAspherJsonFile = false;



            if (find(options.begin(), options.end(), "-f") != options.end())
                shelxFreeFormat = true;
            if (find(options.begin(), options.end(), "-e") != options.end())
                electronScatteringIfNoAspherJsonFile = true;
            if (find(options.begin(), options.end(), "-o") != options.end())
                newImplementation = false;


            bool unitCellChargeIfNoAspherJsonFile = false;
            double unitCellCharge = 0;
            if (optionsWithValues.find("-charge") != optionsWithValues.end())
            {
                unitCellChargeIfNoAspherJsonFile = true;
                unitCellCharge = stod(optionsWithValues["-charge"]);
            }


            if (!filesystem::exists(structureFilePath))
                on_error::throwException("declared structure file missing in the current direcory", __FILE__, __LINE__);

            string suffix = structureFilePath.extension().string();

            int nTwinComponents = 1;
            Matrix3i twinMatrix;

            

            read_structure(structureFile, crystal, nTwinComponents, twinMatrix);


            int nSymm = crystal.spaceGroup.nSymmetryOperationsInSubset();
            vector<Matrix3i> rotations(nSymm);

            for (int i = 0; i < nSymm; i++)
                crystal.spaceGroup.getSpaceGroupOperation(0, 0, i).getRotation(rotations[i]);
            for (int i = 0; i < nSymm; i++)
                rotations.push_back(-1 * rotations[i]);


            vector<Vector3i> hkls, hklAll;

            readHkl(hklFile, hkls, shelxFreeFormat);

            makeWholeHklSet(hkls, rotations, twinMatrix, nTwinComponents, hklAll);

            WallClockTimer timer;

            timer.start();
            // make form factors calculator
            header += string("Uses files ") + structureFile + string(" & ") + hklFile + string("\n\n");

            clog << header
                << "\n" << "file created at " << ctime(&time_now) << "\n";

            

            
            bool aspherJsonPresent;
            auto calculator = 
                sfCalculatorFromJsonFile(crystal, aspherJsonPresent, electronScatteringIfNoAspherJsonFile, 
                                         unitCellChargeIfNoAspherJsonFile, unitCellCharge, electronScattering, jobName);
            

            // is it taam2tsc? i.e. discamb2tsc run without arguments and making default TAAM calculations

            string modelName = calculator->name();
            if (modelName == string("MATTS") && cifAndHklFromDirectorySearch)
                defaultTaamRun = true;

            //-----------
            // generate the output file



            //string stem = structureFilePath.stem().string();
            string outputFileName = jobName + string(".tsc");
            ofstream out(outputFileName);


            
            out << "TITLE: "; 
            if (electronScattering)
                out << "electron ";
            else
                out << "X-ray ";

            out << "form factors for " << jobName << " generated with discamb2tsc\n";

            out << "AD: FALSE\nSYMM: expanded\nSCATTERERS:";
            for (auto const& atom : crystal.atoms)
                out << " " << atom.label;
            out << "\n";
            vector<pair<string, string> > modelInfo;

            calculator->getModelInformation(modelInfo);
            out << "FORM FACTORS SOURCE:\n"
                << "    SOFTWARE - discambMATTS2tsc\n";
                
            for (auto& item : modelInfo)
                out << "    " << item.first << " - " << item.second << "\n";

            out << "DATA:\n";

            timer.start();

            vector< vector<complex<double> > > formFactors;
            nAtoms = crystal.atoms.size();
            vector<bool> includeAtom(nAtoms, true);
            calculator->calculateFormFactors(hklAll, formFactors, includeAtom);
            if (!defaultTaamRun)
                clog << "calculation of form factors: " << timer.stop() << " ms\n";
            for (int hklIdx = 0; hklIdx < hklAll.size(); hklIdx++)
            {
                out << hklAll[hklIdx][0] << " " << hklAll[hklIdx][1] << " " << hklAll[hklIdx][2];
                for (auto const& f : formFactors[hklIdx])
                    out << " " << setprecision(6) << fixed << f.real() << "," << f.imag();
                out << "\n";
            }



            out.close();

            //################### discamb2tsc.log ################



            //################ end of discamb2tsc.log ############

            addMattsCitation(discambMATTS2tscVersion);

            std::clog.rdbuf(clog_orginal_rdbuf);
            logFile.close();

        }
        catch (nlohmann::json::parse_error& e)
        {
            std::cout << "message: " << e.what() << '\n'
                << "exception id: " << e.id << '\n'
                << "byte position of error: " << e.byte << std::endl;

            stringstream ss;

            ss << "message: " << e.what() << '\n'
               << "exception id: " << e.id << '\n'
               << "byte position of error: " << e.byte << std::endl;

            string message = "Error when parsing JSON file\n" + ss.str();
            on_error::throwException(message, __FILE__, __LINE__);
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
        ofstream out("discamb_error.log" , ostream::app);
        out << e.what();
        out.close();
    }
}


