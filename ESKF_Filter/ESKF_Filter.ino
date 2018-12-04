#include "LSM6DSM.h"
#include <Wire.h>
#include <math.h>

//LSM6DSM definitions
#define INTERRUPT_PIN 2  // interrupt1 pin definitions, data ready

/* Specify sensor parameters (sample rate is twice the bandwidth)
 * choices are:
 AFS_2G, AFS_4G, AFS_8G, AFS_16G  
 GFS_245DPS, GFS_500DPS, GFS_1000DPS, GFS_2000DPS 
 AODR_12_5Hz, AODR_26Hz, AODR_52Hz, AODR_104Hz, AODR_208Hz, AODR_416Hz, AODR_833Hz, AODR_1660Hz, AODR_3330Hz, AODR_6660Hz
 GODR_12_5Hz, GODR_26Hz, GODR_52Hz, GODR_104Hz, GODR_208Hz, GODR_416Hz, GODR_833Hz, GODR_1660Hz, GODR_3330Hz, GODR_6660Hz
 */ 
static LSM6DSM::Ascale_t Ascale = LSM6DSM::AFS_4G;
static LSM6DSM::Gscale_t Gscale = LSM6DSM::GFS_245DPS;
static LSM6DSM::Rate_t AODR = LSM6DSM::ODR_1660Hz;
static LSM6DSM::Rate_t GODR = LSM6DSM::ODR_1660Hz;

// Gyro and accel parameters
static float GYRO_NOISE_DEN  = 3.8; // [mdps/sqrt(Hz)]
static float ACCEL_NOISE_DEN = 80;  // [ug/sqrt(Hz)] (Different value depending on the FS)

// Random Bias to initiate the object
float ACCEL_BIAS[3] = {0.0, 0.0, 0.0};
float GYRO_BIAS[3]  = {0.0, 0.0, 0.0};

static LSM6DSM lsm6dsm(Ascale, Gscale, AODR, GODR, ACCEL_BIAS, GYRO_BIAS);

//Declaring filter variables
// States
float q[4] = {1.0};  // Quaternion initialization

// Declare error state
float dq[3]; // Error-State vector
float P[9];  // Covariance square matrix. Order, from top to bottom, left to right.

uint32_t t_last = micros();

void setup() 
{
    Serial.begin(115200);
   
    // Configure interrupt
    pinMode(INTERRUPT_PIN, INPUT);

    // Start I^2C 
    Wire.begin(TWI_PINS_20_21);           // set master mode 
    Wire.setClock(400000);  // I2C frequency at 400 kHz  
    delay(100);
    
    lsm6dsm.begin();

    calibrate();
    
    // Gyro std. dev: ("Rate noise density")*sqrt("Sampling Rate"/2) [rad/s];
    float gyro_stdDev =  sqrtf(1660/2)*GYRO_NOISE_DEN; // [mdeg/s]
    gyro_stdDev = gyro_stdDev*M_PI/180000;
    float gyro_noise = gyro_stdDev*gyro_stdDev;
    
    float accel_stdDev =  sqrtf(1660/2)*ACCEL_NOISE_DEN; // [ug/s]
    accel_stdDev = accel_stdDev/1000000; // [g/s]
    float accel_noise = accel_stdDev*accel_stdDev;

              float b[4] = {1.0, -4.0, 3.0, 2.0};
    normalizeVector(4,b,b);
    for (int ii=0; ii<4; ++ii)
    {
      Serial.println(b[ii],4);
    }
    delay(5000);

}

void loop() 
{
  float Q;
  if (lsm6dsm.checkNewData()) 
  {  
    float ax=0, ay=0, az=0, gx=0, gy=0, gz=0;

    // acc [g], gyro [deg/s]
    lsm6dsm.readData(ax, ay, az, gx, gy, gz);
    
    uint32_t t_delta = t_last - micros();
    t_last = micros();
    
    // Update state estimation
    float Fq[16];
    computeFq(&gx, &gy, &gz, t_delta, Fq);
    matrixByVector(4,4,Fq,q,q);
    normalizeVector(4,q,q);

    // Update error-state estimation
    float Fdq[3]; 

  }
}

void calibrate()
{
  Serial.println("DO NOT move the IMU. Calibration starting in:");
  Serial.println("3");
  delay(1000);
  Serial.println("2");
  delay(1000);
  Serial.println("1");
  delay(1000);
  Serial.println("Calibrating...");
  
  lsm6dsm.calibrate(GYRO_BIAS, ACCEL_BIAS);

  Serial.print("Calibrated");
}

void computeFq(float *gx, float *gy, float *gz, uint32_t t_delta, float *Fx)
{
  // First Column
  Fx[0]  = 1.0;
  Fx[1]  = (*gx)*t_delta/2.0;
  Fx[2]  = (*gy)*t_delta/2.0;
  Fx[3]  = (*gz)*t_delta/2.0;
  // Second Column
  Fx[4]  = -(*gx)*t_delta/2.0;
  Fx[5]  =  1.0;
  Fx[6]  = -(*gz)*t_delta/2.0;
  Fx[7]  =  (*gy)*t_delta/2.0;
  // Third Column
  Fx[8]  = -(*gy)*t_delta/2.0;
  Fx[9]  =  (*gz)*t_delta/2.0;
  Fx[10] =  1.0;
  Fx[11] = -(*gx)*t_delta/2.0;
  // Fourth Column
  Fx[12] = -(*gz)*t_delta/2.0;
  Fx[13] = -(*gy)*t_delta/2.0;
  Fx[14] =  (*gx)*t_delta/2.0;
  Fx[15] =  1.0;
}

void matrixByVector(uint8_t nrow, uint8_t ncol, float *mat, float *vec, float *res)
{
  float res_t[ncol];
  for (int ii = 0; ii < nrow; ++ii)
  {
    for (int jj = 0; jj < ncol; ++jj)
    {
      res_t[ii] += mat[ii+nrow*jj]*vec[jj];
    }
  }

  for (int jj = 0; jj < ncol; ++jj)
  {
    res[jj] = res_t[jj];
  }
}

void normalizeVector(int nele, float* vec, float *res)
{
  float norm;
  for (int ii=0; ii<nele; ++ii)
  {
    norm += vec[ii]*vec[ii];
  }
  
  for (int ii=0; ii<nele; ++ii)
  {
    res[ii] = vec [ii]/sqrtf(norm);
  }
  
}

