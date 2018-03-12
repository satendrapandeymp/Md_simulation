using namespace std;
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>
#include <fstream>
#include <stdlib.h>
#include "Declearation.h"
#include "Multiparticle.h"

/////////////////////////////////////////////////////////////////////////////################################################################///////////////////////////////////////////////////////////

// Writing functions here

// Update distances
void distance_cal(float pos[864][3], float * dist);

// To calculate force Brute_force for now
float force_cal(float posx, float posy, float posz, float * accx, float * accy, float * accz);

// To update pos and velocity
float update_all(float * posx, float * posy, float * posz, float * velx, float  * vely, float  * velz, float accx, float accy, float accz,float last_accx, float last_accy, float last_accz, float delta_time);

// For periodic Boundry Condition
void periodic_bc();

// To calculate total energy
void potential_energy ();
void kinetic_energy ();

// To calculate pair dist functions
void pair_dist();

// To calculate velocity dist function
void vel_auto();

// clearing matrix
void clearing();

/////////////////////////////////////////////////////////////////////////////################################################################///////////////////////////////////////////////////////////

int main()
{
    float potential, kinetic;

    // Reading from file
    declearation("data.txt");

    // calculating dist for further use
    //distance_cal(pos, dist);

    // To calculate force on each atom
    for ( i = 0; i<864; i++)
    {
      potential = force_cal(pos[i][0], pos[i][1], pos[i][2], &acc[i][0], &acc[i][1], &acc[i][2]);
    }

    // Now starting the real loop
    // for 100 time steps
    for (int step = 0; step < 100; step++)
    {
      potential = 0;
      kinetic = 0;

      // In each step going in loop to update everything
      for (i = 0; i<864; i++)
      {
            last_acc[i][0] =acc[i][0];
            last_acc[i][1] =acc[i][1];
            last_acc[i][2] =acc[i][2];

            // To calculate force on each atom
            potential = potential + force_cal(pos[i][0], pos[i][1], pos[i][2], &acc[i][0], &acc[i][1], &acc[i][2]);

      }

      for (i = 0; i<864; i++)
      {

            // Updating position and velocity
            kinetic = kinetic + update_all(&pos[i][0], &pos[i][1], &pos[i][2], &vel[i][0], &vel[i][1] ,  &vel[i][2] , acc[i][0], acc[i][1], acc[i][2],last_acc[i][0], last_acc[i][1], last_acc[i][2], delta_time);
      }

      for (i = 0; i<1; i++)
      {
            cout << pos[i][0] << "  " << pos[i][1] << "  " << pos[i][2] << "  " << (kinetic + potential) << endl;
      }

      clearing();

    }




return 0;

}

////////////////////////////////////////////////////////////////////////####################################################################///////////////////////////////////////////////////////////

// Will write functions here


float update_all(float * posx, float * posy, float * posz, float * velx, float  * vely, float  * velz, float accx, float accy, float accz,float last_accx, float last_accy, float last_accz, float delta_time)
{
  float kinetic;
  float mbk = 1/.0247713;

  *posx = *posx + *velx * delta_time + accx * delta_time * delta_time/2;
  *posy = *posy + *vely * delta_time + accy * delta_time * delta_time/2;
  *posz = *posz + *velz * delta_time + accz * delta_time * delta_time/2;

  // update position with periodic box condition
  if (*posx > 34.78) *posx = *posx - 34.78;
  if (*posx < 0) *posx = *posx + 34.78;
  if (*posy > 34.78) *posy = *posy - 34.78;
  if (*posy < 0) *posy = *posy + 34.78;
  if (*posz > 34.78) *posz = *posz - 34.78;
  if (*posz < 0) *posz = *posz + 34.78;

  kinetic = mbk * (pow(*velx, 2) + pow(*vely, 2) + pow(*velz, 2))/2;

  *velx = *velx + (accx + last_accx) * delta_time/2;
  *vely = *vely + (accy + last_accy) * delta_time/2;
  *velz = *velz + (accz + last_accz) * delta_time/2;

  return kinetic;

}

// clearing attraction matrix
void clearing()
{
  for (int i=0; i<864; i++)
  {
    for (int j = 0; j<3; j++)
    {
      acc[i][j] = 0;
    }
  }
}







// Don't mess with what's written below :p
//*********************************************************************************************************************************************************
//---------------------------------------------------------------------------------------------------------------------------------------------------------
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Update pos
void distance_cal(float pos[864][3], float * dist)
{
  for (int i = 0; i<864; i++)
  {
    float temp =  pow( (pow(pos[i][0],2) + pow(pos[i][1],2) + pow(pos[i][2],2) ),.5);
    dist[i] = temp;
  }
}



// Force calculation
float force_cal(float posx, float posy, float posz, float * accx, float * accy, float * accz)
{
  float potential = 0;
// For all Atoms
  int flag = 0;
  if ((posx-17.39) > 9.8) flag++;
  if ((posy-17.39) > 9.8) flag++;
  if ((posz-17.39) > 9.8) flag++;

for (int i = 0; i<27; i++)
{

  // Mapping to 3D
  int num = i;
  int cor_x, cor_y, cor_z;
  cor_x = (num/9)-1;
  cor_y = ((num%9)/3)-1;
  cor_z = (num%3)-1;

  //cout << cor_x << cor_y << cor_z << endl;

if ((abs(cor_x) + abs(cor_y) + abs(cor_z)) <= flag)
{

  for (int j = 0; j<864; j++)
          {
              float r = pow( (pow( (pos[j][0] + cor_x * 34.78 )- posx, 2) + pow((pos[j][1] + cor_y * 34.78 )- posy, 2) + pow((pos[j][2] + cor_z * 34.78 )- posz, 2) ), .5);
              if (r < 7.65 && r > .1)
              {
                  // Beware of units, its acc not force
                  // letus measure energy in form of kt
                  potential = potential + 4 * 120 * ((pow((sigma/r),12))  - (pow((sigma/r),6)));
                  float force = 4 * kbm * ((pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                  *accx = *accx + force * ((pos[j][0] + cor_x * 34.78)-posx)/r;
                  *accy = *accy + force * ((pos[j][1] + cor_x * 34.78)-posy)/r;
                  *accz = *accz + force * ((pos[j][2] + cor_x * 34.78 )-posz)/r;

              }
            }

}

}

return potential;


  // For stoms at sides periodic boundry condition Too messy and can't do in this much short time
  //////////////////////////////////////##################################For one edge//////////////////////////////////////***********************************//////////////////////////////////////

/*
  if (flag == 1)
        {
          if (abs(posx - 17.389) > 10 )
          {


              if ((posx - 17.389) > 10 ))
              {


                for (int i = 0; i<864; i++)
                        {
                          if ((pos[i][0] +  34.78) - posx) < 7.65)
                          {
                            float r = pow( (pow((pos[i][0] + 34.78) - posx, 2) + pow(pos[i][1] - posy, 2) + pow(pos[i][2]-posz, 2) ), .5);
                            if (r < 7.65 && r > .1)
                            {
                                // Beware of units, its acc not force
                                float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                *accx = *accx + force * ((pos[i][0] + 34.78) - posx);
                                *accy = *accy + force * (pos[i][1]-posy);
                                *accz = *accz + force * (pos[i][2]-posz);

                                //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                              }
                            }
                          }

              }
              else
              {

                for (int i = 0; i<864; i++)
                        {
                          if ((pos[i][0] - posx) > 27)
                          {
                            float r = pow( (pow((pos[i][0] - 34.78) - posx, 2) + pow(pos[i][1] - posy, 2) + pow(pos[i][2]-posz, 2) ), .5);
                            if (r < 7.65 && r > .1)
                            {
                                // Beware of units, its acc not force
                                float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                *accx = *accx + force * ((pos[i][0] - 34.78) - posx);
                                *accy = *accy + force * (pos[i][1]-posy);
                                *accz = *accz + force * (pos[i][2]-posz);

                                //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                              }
                            }
                          }

              }
          }

          // For y side-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          if (abs(posy - 17.389) > 10 )
          {


              if ((17.389 - posy) > 10 ))
              {

                for (int i = 0; i<864; i++)
                        {
                          if ((pos[i][1] - posy) > 27)
                          {
                            float r = pow( (pow((pos[i][0]) - posx, 2) + pow((pos[i][1]- 34.78) - posy, 2) + pow(pos[i][2]-posz, 2) ), .5);
                            if (r < 7.65 && r > .1)
                            {
                                // Beware of units, its acc not force
                                float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                *accx = *accx + force * (pos[i][0] - posx);
                                *accy = *accy + force * ((pos[i][1] - 34.78) - posy);
                                *accz = *accz + force * (pos[i][2]-posz);

                                //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                              }
                            }
                          }

              }
              else
              {

                for (int i = 0; i<864; i++)
                        {
                          if ((posy - pos[i][1]) > 27)
                          {
                            float r = pow( (pow((pos[i][0]) - posx, 2) + pow((pos[i][1] + 34.78) - posy, 2) + pow(pos[i][2]-posz, 2) ), .5);
                            if (r < 7.65 && r > .1)
                            {
                                // Beware of units, its acc not force
                                float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                *accx = *accx + force * (pos[i][0] - posx);
                                *accy = *accy + force * ((pos[i][1] + 34.78) - posy);
                                *accz = *accz + force * (pos[i][2]-posz);

                                //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                              }
                            }
                          }

              }
          }

          // for z ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


          if (abs(posz - 17.389) > 10 )


                        if ((17.389 -  posz ) > 10 ))
                        {


                          for (int i = 0; i<864; i++)
                                  {
                                    if ((pos[i][2] - posz) > 27)
                                    {
                                      float r = pow( (pow((pos[i][0]) - posx, 2) + pow((pos[i][1]) - posy, 2) + pow((pos[i][2]- 34.78)-posz, 2) ), .5);
                                      if (r < 7.65 && r > .1)
                                      {
                                          // Beware of units, its acc not force
                                          float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                          *accx = *accx + force * (pos[i][0] - posx);
                                          *accy = *accy + force * (pos[i][1] - posy);
                                          *accz = *accz + force * ((pos[i][2] - 34.78) - posz);

                                          //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                          //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                          //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                                        }
                                      }
                                    }

                        }
                        else
                        {

                          for (int i = 0; i<864; i++)
                                  {
                                    if ((posy - pos[i][1]) < 27)
                                    {
                                      float r = pow( (pow((pos[i][0]) - posx, 2) + pow((pos[i][1]) - posy, 2) + pow((pos[i][2]- 34.78)-posz, 2) ), .5);
                                      if (r < 7.65 && r > .1)
                                      {
                                          // Beware of units, its acc not force
                                          float force = 4 * kbm * ( (pow((sigma/r),13) * 12)  - (pow((sigma/r),7) * 6))/ sigma;
                                          *accx = *accx + force * (pos[i][0] - posx);
                                          *accy = *accy + force * ((pos[i][1] - 34.78) - posy);
                                          *accz = *accz + force * (pos[i][2]-posz);

                                          //float check = pow( pow((pos[i][0]-posx)/r, 2) + pow((pos[i][1]-posx)/r, 2) + pow((pos[i][2]-posx)/r, 2), .5);
                                          //cout << check << " " << r << " " << posx <<" "  << posy << " "  << posz <<" "  << pos[i][0] << " "  << pos[i][1] <<" "  << pos[i][2] << endl;

                                          //cout << pow(*accx, 2) + pow(*accy, 2) + pow(*accz, 2)  << "  " <<  force * force << endl;
                                        }
                                      }
                                    }

                        }

        }

        */

}
