#include <fstream>
#include <iostream>
#include <stdio.h>
#include "../include/mtpParam.h"


namespace matersdk {
namespace mtpr {

MtpParam::MtpParam() {}

MtpParam::MtpParam(const std::string& filename)
{
    this->_load(filename);
}

void MtpParam::_load(const std::string& filename)
{
    this->_alpha_count = 0;
    
    std::ifstream ifs(filename);
    if ( !ifs.is_open() ) 
        MtpError((std::string)"Cannot open " + filename);

    char tmpline[1000];
    std::string tmpstr;

    ifs.getline(tmpline, 1000);
    int len = (int)((std::string)tmpline).length();
    if (tmpline[len-1] == '\r') // Ensures compatibility between Linux and Windows line endings
        tmpline[len-1] = '\0';
    if ((std::string)tmpline != "MTP") 
        MtpError("Can read only MTP potentials");


    // version reading block
    ifs.getline(tmpline, 1000);
    len = (int)((std::string)tmpline).length();
    if (tmpline[len-1] == '\r')
        tmpline[len-1] = '\0';
    if (((std::string)tmpline) != "version = 1.1.0")
        MtpError("MTP file must have version 1.1.0");

    // name/description reading block
    ifs >> tmpstr; 
    if (tmpstr == "potential_name") // optional
    {
        ifs.ignore(2);
        ifs >> this->pot_desc;
        ifs >> tmpstr;
    }

    if (tmpstr == "scaling")    // optional
    {
        ifs.ignore(2);
        ifs >> this->scaling;
        ifs >> tmpstr;
    }

    if (tmpstr == "species_count")
    {
        ifs.ignore(2);
        ifs >> this->species_count;
        ifs >> tmpstr;
    }
    else
        species_count = 0;
    
    if (tmpstr == "potential_tag")
    {
        std::getline(ifs, tmpstr);
        ifs >> tmpstr;
    }

    if (tmpstr != "radial_basis_type")
        MtpError("Error reading .mtp file " + filename);
    ifs.ignore(2);
    ifs >> this->rbasis_type;

    if (rbasis_type == "BChebyshev") {
        //p_RadialBasis = new Basis_Chebyshev(ifs);
    }
    else if (rbasis_type == "BChebyshev_repuls") {
        //p_RadialBasis = new Basis_Chebyshev_repuls(ifs);
    }
    else if (rbasis_type == "RBChebyshev") {
        //p_RadialBasis = new RB_Chebyshev(ifs);
    }
    else if (rbasis_type == "RBChebyshev_repuls") {
        //p_RadialBasis = new RB_Chebyshev_repuls(ifs);
    }
    else if (rbasis_type == "RBShapeev") {
        //p_RadialBasis = new Basis_Shapeev(ifs);
        //Warning("Non-linear MTP is not fully compatible with RBShapeev");
    }
    else if (rbasis_type == "RBTaylor") {
        //p_RadialBasis = new Basis_Taylor(ifs);
    }
    else
        MtpError("Wrong radial basis type");

    // We do not need double scaling
    //if (p_RadialBasis->scaling != 1.0) {
    //    scaling *= p_RadialBasis->scaling;
    //    p_RadialBasis->scaling = 1.0;
    //}

    ifs >> tmpstr;
    if (tmpstr == "magnetic_basis_type") {
        ifs.ignore(2);
        ifs >> this->mbasis_type;
        if (mbasis_type == "BChebyshev") {
            //p_MagneticBasis = new Basis_Chebyshev(ifs);
        }
        else if (mbasis_type == "BChebyshev_repuls") {
            //p_MagneticBasis = new Basis_Chebyshev_repuls(ifs);
        }
        else if (mbasis_type == "RBChebyshev") {
            //p_MagneticBasis = new RB_Chebyshev(ifs);
        }
        else if (mbasis_type == "RBChebyshev_repuls") {
            //p_MagneticBasis = new RB_Chebyshev_repuls(ifs);
        }
        else if (mbasis_type == "RBShapeev") {
            //p_MagneticBasis = new Basis_Shapeev(ifs);
        }
        else if (mbasis_type == "RBTaylor") {
            //p_MagneticBasis = new Basis_Taylor(ifs);
        }
        else if (mbasis_type == "") {
            ;
        }
        else
            MtpError("Wrong magnetic basis type");

        //mbasis_size = p_MagneticBasis->size;

        // We do not need double scaling
        /*if (p_MagneticBasis->scaling != 1.0) {
            scaling *= p_MagneticBasis->scaling;
            p_MagneticBasis->scaling = 1.0;
        }*/
        ifs >> tmpstr;
    }

    ifs.getline(tmpline, 1000);  // min_dist
    ifs.getline(tmpline, 1000);  // max_dist
    ifs.getline(tmpline, 1000);  // radial_basis_size
    ifs >> tmpstr; // radial_funcs_count
    if (tmpstr != "radial_funcs_count") 
        MtpError("Error reading .mtp file");
    ifs.ignore(2);
    ifs >> this->radial_func_count;

    // Radial coeffs initialization
    //int pairs_count = species_count * species_count;
    //int mbasis_size2 = mbasis_size * mbasis_size;
    
    ifs >> tmpstr;
    if (tmpstr != "alpha_moments_count")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(2);
    ifs >> this->_alpha_moments_count;
    if (ifs.fail())
        MtpError("Error reading .mtp file "+filename);

    ifs >> tmpstr;
    if (tmpstr != "alpha_index_basic_count")
        MtpError("Error reading .mtp file.");
    ifs.ignore(2);
    ifs >> this->_alpha_index_basic_count;
    if (ifs.fail())
        MtpError("Error reading .mtp file "+filename);
    
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_basic")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(4);
    this->_alpha_index_basic = new int[this->_alpha_index_basic_count][4];
    if (this->_alpha_index_basic == nullptr)
        MtpError("Memory allocation error");
    
    int radial_func_max = -1;
    for (int i=0; i<this->_alpha_index_basic_count; i++) 
    {
        char tmpch;
        ifs.ignore(1000, '{');
        ifs >> this->_alpha_index_basic[i][0] >> tmpch 
            >> this->_alpha_index_basic[i][1] >> tmpch 
            >> this->_alpha_index_basic[i][2] >> tmpch
            >> this->_alpha_index_basic[i][3];
        if (ifs.fail())
            MtpError("Error reading .mtp file "+filename);
        
        if (this->_alpha_index_basic[i][0] > radial_func_max)
            radial_func_max = this->_alpha_index_basic[i][0];
    }

    if (radial_func_max != radial_func_count-1)
        MtpError("Wrong number of radial functions specified in "+filename);

    ifs.ignore(1000, '\n');
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_times_count")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(2);
    ifs >> this->_alpha_index_times_count;
    if (ifs.fail())
        MtpError("Error reading .mtp file "+filename);
    
    ifs >> tmpstr;
    if (tmpstr != "alpha_index_times")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(4);

    this->_alpha_index_times = new int[this->_alpha_index_times_count][4];
    if (this->_alpha_index_times == nullptr)
        MtpError("Memory allocation error "+filename);
    
    for (int i=0; i<this->_alpha_index_times_count; i++)
    {
        char tmpch;
        ifs.ignore(1000, '{');
        ifs >> this->_alpha_index_times[i][0] >> tmpch
            >> this->_alpha_index_times[i][1] >> tmpch
            >> this->_alpha_index_times[i][2] >> tmpch
            >> this->_alpha_index_times[i][3];
        if (ifs.fail())
            MtpError("Error reading .mtp file "+filename);
    }

    ifs.ignore(1000, '\n');

    ifs >> tmpstr;
    if (tmpstr != "alpha_scalar_moments")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(2);
    ifs >> this->_alpha_scalar_moments;
    if (this->_alpha_scalar_moments < 0)
        MtpError("Error reading .mtp file "+filename);
    
    this->_alpha_moment_mapping = new int[this->_alpha_scalar_moments];
    if (this->_alpha_moment_mapping == nullptr)
        MtpError("Memory allocation error");

    ifs >> tmpstr;
    if (tmpstr != "alpha_moment_mapping")
        MtpError("Error reading .mtp file "+filename);
    ifs.ignore(4);
    for (int i=0; i<this->_alpha_scalar_moments; i++)
    {
        char tmpch = ' ';
        ifs >> this->_alpha_moment_mapping[i] >> tmpch;
        if (ifs.fail())
            MtpError("Error reading .mtp file "+filename);
    }
    ifs.ignore(1000, '\n');
}

MtpParam::~MtpParam() {
    delete this->_alpha_index_basic;
    delete this->_alpha_index_times;
    delete this->_alpha_moment_mapping;
}

};  // namespace : mtpr
};  // namespace : matersdk