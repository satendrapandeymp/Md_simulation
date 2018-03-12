#ifndef MULTIPARTICLE_H_INCLUDED
#define MULTIPARTICLE_H_INCLUDED


// Declearing 2D array
// Will modify this section
void declearation(string fileName)
{

    // Testing
    string temp;

    ifstream file(fileName.c_str());

    if(file.good())
        {

          for (i = 0; i<864; i++)
            {
                // Read from file for this section
                file >> temp;
                pos[i][0] =(double)::atof(temp.c_str());
                file >> temp;
                pos[i][1] = (double)::atof(temp.c_str());
                file >> temp;
                pos[i][2] = (double)::atof(temp.c_str());

                // This section could be replace by zero
                file >> temp;
                vel[i][0] = (double)::atof(temp.c_str());
                file >> temp;
                vel[i][1] = (double)::atof(temp.c_str());
                file >> temp;
                vel[i][2] = (double)::atof(temp.c_str());
                acc[i][0] = 0;
                acc[i][1] = 0;
                acc[i][2] = 0;

            }
        }

    else
        {
            std::cout << "Can't find the file" << endl;
        }

}


#endif
