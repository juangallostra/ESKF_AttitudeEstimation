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
static LSM6DSM::Gscale_t Gscale = LSM6DSM::GFS_2000DPS;
static LSM6DSM::Rate_t AODR = LSM6DSM::ODR_1660Hz;
static LSM6DSM::Rate_t GODR = LSM6DSM::ODR_1660Hz;

// Gyro and accel parameters
static float GYRO_NOISE_DEN  = 3.8f; // [mdps/sqrt(Hz)]
static float ACCEL_NOISE_DEN = 80.0f;  // [ug/sqrt(Hz)] (Different value depending on the FS)

// Random Bias to initiate the object
float ACCEL_BIAS[3] = {0.0f, 0.0f, 0.0f};
float GYRO_BIAS[3]  = {0.0f, 0.0f, 0.0f};

static LSM6DSM lsm6dsm(Ascale, Gscale, AODR, GODR, ACCEL_BIAS, GYRO_BIAS);

//Declaring filter variables
// States
float q_[4] = {1.0f,0.0f,0.0f,0.0f};  // Quaternion initialization
float dq_[3] = {0.0f,0.0f,0.0f}; // Error-State vector
float P_[9] = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};  // Covariance square matrix. Order, from top to bottom, left to right.
float gyro_noise, accel_noise;
uint32_t t_last = micros();

uint16_t loop_counter = 0;
    
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

    delay(1000);
    
    // Gyro std. dev: ("Rate noise density")*sqrt("Sampling Rate"/2) [rad/s];
    float gyro_stdDev =  sqrtf(1660/2)*GYRO_NOISE_DEN; // [mdeg/s]
    gyro_stdDev = gyro_stdDev*M_PI/180000; // [rad/s]
    // gyro_stdDev = 0.0024;
    gyro_noise = gyro_stdDev*gyro_stdDev;
    
    
    float accel_stdDev =  sqrtf(1660/2)*ACCEL_NOISE_DEN; // [ug/s]
    accel_stdDev = accel_stdDev/1000000; // [g/s]
    accel_stdDev = 0.10283;
    accel_noise = accel_stdDev*accel_stdDev;
    
    // Initialize Covariance
    P_[0] = 0.0000000001f;
    P_[1] = 0.0f;
    P_[2] = 0.0f;
    P_[3] = 0.0f;
    P_[4] = 0.0000000001f;
    P_[5] = 0.0f;
    P_[6] = 0.0f;
    P_[7] = 0.0f;
    P_[8] = 0.0000000001f;

    q_[0] = 1.0f;
    q_[1] = 0.0f;
    q_[2] = 0.0f;
    q_[3] = 0.0f;

    dq_[0] = 0.0;
    dq_[1] = 0.0;
    dq_[2] = 0.0;

   /* Serial.println("Gyro Noise");
    Serial.println(gyro_noise,10);
    Serial.println("Accel Noise");
    Serial.println(accel_noise,10);*/
}

void loop() 
{
  float Fq[16] = {0.0f};
  float Fdq[9] = {0.0f};
  float Q[9]   = {0.0f};
  float Hdq[9] = {0.0f};
  float Na[9]  = {0.0f};  
  float Z[9]   = {0.0f};  
  float Zinv[9]= {0.0f};  
  float K[9] = {0.0f};
  float Hdq_tr[9]= {0.0f};  
  float R[9]= {0.0f};  
  float R_tr[9]= {0.0f};  
  float acc_est[3]= {0.0f};  
  float z[3]= {0.0f};  
  float aux[9]= {0.0f};  
  float I3[9] = {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f};
  float S[9]= {0.0f};
  float qL[16]= {0.0f};
  float dq_aux[4] = {1.0f};
  float G[9]= {0.0f};  
  float euler[3]= {0.0f};

  
//if (loop_counter <=255)
//{
  //++loop_counter;
  //printMatrix3("K at the beggining",K);
  if (lsm6dsm.checkNewData()) 
  {  
    /*Serial.println(" ");
    Serial.println(" ");
    Serial.println(" ");
    Serial.println(" ");
    Serial.print("LOOP COUNTER:   ");
    Serial.println(loop_counter);*/
        
    float ax=0.0f, ay=0.0f, az=0.0f, gx=0.0f, gy=0.0f, gz=0.0f;
    float gravity[3] = {0.0f, 0.0f, 9.80665f};
    
    // acc [g], gyro [deg/s]
    lsm6dsm.readData(ax, ay, az, gx, gy, gz);
    gx = gx*M_PI/180;
    gy = gy*M_PI/180;
    gz = gz*M_PI/180;
    
    ax = ax*gravity[2];
    ay = ay*gravity[2];
    az = az*gravity[2];

    ax = int(ax*10000)/10000.0  ;
    ay = int(ay*10000)/10000.0  ;
    az = int(az*10000)/10000.0  ;
    gx = int(gx*10000)/10000.0  ;
    gy = int(gy*10000)/10000.0  ;
    gz = int(gz*10000)/10000.0  ;
    
    float gyro_m[3] = {gx, gy, gz};
    float acc_m[3] = {ax, ay, az};

    //printVector3("Gyrometer measure", gyro_m);
   // printVector3("Accel measure", acc_m);
            
    uint32_t t_delta = micros() - t_last;
    t_last = micros();
    
    float t_deltaf = t_delta / 1000000.0;
    /*Serial.println("Delta time:");
    Serial.println(t_deltaf, 12);*/
        
    // Update state estimation
    computeFq(&gx, &gy, &gz, t_deltaf, Fq);
    matrixByVector4(Fq,q_,q_);
    normalizeVector4(q_,q_);

   // printMatrix4("Jacobian Fq:", Fq);
    //printVector4("Quaternion estimation update:", q_);
  
    // Update error-state estimation
    computeFdq(&gx, &gy, &gz, t_deltaf, Fdq);
    computeQ(t_deltaf, Q);
    matrixByVector3(Fdq,dq_, dq_);
    propagateCovariance(P_,Fdq,P_);
    addMatrices(3,3,P_,Q,P_);

    //printMatrix3("Jacobian Fdq", Fdq);
    //printMatrix3("Covariance P", P_);
    ///printMatrix3("Noise matrix Q", Q);
    //printVector3("Error state dq", dq_);      
    
    // Correction
    computeHdq(q_,Hdq);  
    computeNa(Na);
    propagateCovariance(P_,Hdq,Z);
    addMatrices(3,3,Z,Na,Z);
    inverseMatrix3(Z, Zinv);
    transpose3(Hdq, Hdq_tr);
    matrixByMatrix(3,P_,Hdq_tr,K);
    matrixByMatrix(3,K,Zinv,K);
      
    /*printMatrix3("Jacobian Hdq", Hdq);
    printMatrix3("Matrix Z", Z);
    printMatrix3("Matrix Zinv", Zinv);
    printMatrix3("Gain K", K);
    printMatrix3("Noise Na", Na);*/
        
    fromqtoR(q_,R);
    transpose3(R,R_tr);
    matrixByVector3(R_tr,gravity,acc_est);
    normalizeVector3(acc_est, acc_est);
    normalizeVector3(acc_m, acc_m);
    z[0] = acc_m[0] - acc_est[0];
    z[1] = acc_m[1] - acc_est[1];
    z[2] = acc_m[2] - acc_est[2];
    matrixByVector3(K,z,dq_);
    
    /*printVector3("Estimated acceler:", acc_est);
    printVector3("Error vector z", z);
    printVector3("Updated error state dq", dq_);*/
    
    matrixByMatrix(3, K, Hdq, aux);
    subtractMatrices(3, 3, I3, aux, S);
    propagateCovariance(P_,S,P_);
    propagateCovariance(Z,K, aux);
    addMatrices(3,3,P_,aux,P_);

    //printMatrix3("S matrix", S);
    //printMatrix3("Covariance P", P_);
    
    leftQuaternion(q_, qL);
    dq_aux [0] = 1.0;
    dq_aux [1] = dq_[0]*0.5; 
    dq_aux [2] = dq_[1]*0.5; 
    dq_aux [3] = dq_[2]*0.5;
    
    matrixByVector4(qL, dq_aux, q_);
    dq_[0] = 0.5*dq_[0];
    dq_[1] = 0.5*dq_[1];
    dq_[2] = 0.5*dq_[2];
    skew(dq_, aux);
    
    subtractMatrices(3,3,I3, aux, G);
    propagateCovariance(P_,G,P_);    
    
    dq_[0] = 0.0;
    dq_[1] = 0.0;
    dq_[2] = 0.0;

    truncateNumbers();
    
    //printMatrix3("Last Covariance P", P_);
    //printVector4("Quaternion Last", q_);
    
    computeEulerAngles(q_,euler);
    Serial.print(q_[0]);
    Serial.print(",");
    Serial.print(q_[1]);
    Serial.print(",");
    Serial.print(q_[2]);
    Serial.print(",");
    Serial.println(q_[3]);
  }
//}
}

void truncateNumbers()
{
  P_[0] = int(P_[0]*1000000000000)/1000000000000.0;
  P_[1] = int(P_[1]*1000000000000)/1000000000000.0;
  P_[2] = int(P_[2]*1000000000000)/1000000000000.0;
  P_[3] = int(P_[3]*1000000000000)/1000000000000.0;
  P_[4] = int(P_[4]*1000000000000)/1000000000000.0;
  P_[5] = int(P_[5]*1000000000000)/1000000000000.0;
  P_[6] = int(P_[6]*1000000000000)/1000000000000.0;
  P_[7] = int(P_[7]*1000000000000)/1000000000000.0;
  P_[8] = int(P_[8]*1000000000000)/1000000000000.0;
  
  /*q_[0] = int(q_[0]*100000)/100000.0;
  q_[1] = int(q_[1]*100000)/100000.0;
  q_[2] = int(q_[2]*100000)/100000.0;
  q_[3] = int(q_[3]*100000)/100000.0;*/
  
}

void calibrate()
{
  Serial.println("DO NOT move the IMU. Calibration starting in:");
  Serial.println("3");
  delay(500);
  Serial.println("2");
  delay(500);
  Serial.println("1");
  delay(500);
  
  Serial.println("Calibrating...");
  lsm6dsm.calibrate(GYRO_BIAS, ACCEL_BIAS);
  Serial.println("Calibrated");
}

void computeFq(float *gx, float *gy, float *gz, float t_delta, float Fx[16])
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

void computeFdq(float *gx, float *gy, float *gz, float t_delta, float Fdq[9])
{
  // First Column
  Fdq[0]  =  1.0;
  Fdq[1]  =  (*gz)*t_delta;
  Fdq[2]  = -(*gy)*t_delta;
  // Second Column
  Fdq[3]  = -(*gz)*t_delta;
  Fdq[4]  =  1.0;
  Fdq[5]  =  (*gx)*t_delta;
  // Third Column
  Fdq[6]  =  (*gy)*t_delta;
  Fdq[7]  = -(*gx)*t_delta;
  Fdq[8]  = 1.0;
}


void computeHdq(float q[4], float Hdq[9])
{
  // First Column
  Hdq[0]  =  0.0;
  Hdq[1]  =  q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  Hdq[2]  = - 2*q[0]*q[1] - 2*q[2]*q[3];
  // Second Column
  Hdq[3]  = -q[0]*q[0] + q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  Hdq[4]  =  0.0;
  Hdq[5]  = 2*q[1]*q[3] - 2*q[0]*q[2];
  // Third Column
  Hdq[6]  = 2*q[0]*q[1] + 2*q[2]*q[3];
  Hdq[7]  = 2*q[0]*q[2] - 2*q[1]*q[3];
  Hdq[8]  = 0.0;
}

void computeQ(float t_delta, float Fdq[9])
{
  // First Column
  Fdq[0]  =  gyro_noise*t_delta*t_delta;
  // Second Column
  Fdq[4]  = Fdq[0];
  // Third Column
  Fdq[8]  = Fdq[0];
}

void computeNa(float Na[9])
{
  // First Column
  Na[0]  =  accel_noise;
  // Second Column
  Na[4]  = Na[0];
  // Third Column
  Na[8]  = Na[0];
}

void propagateCovariance(float P[9], float F[9], float res[9])
{
  // Propagates covariance as P = F*P*F'
  float res_t[9];
  // First Column
  res_t[0] =  F[0]*(F[0]*P[0] + F[3]*P[1] + F[6]*P[2]) + F[3]*(F[0]*P[3] + F[3]*P[4] + F[6]*P[5]) + F[6]*(F[0]*P[6] + F[3]*P[7] + F[6]*P[8]);
  res_t[1] =  F[0]*(F[1]*P[0] + F[4]*P[1] + F[7]*P[2]) + F[3]*(F[1]*P[3] + F[4]*P[4] + F[7]*P[5]) + F[6]*(F[1]*P[6] + F[4]*P[7] + F[7]*P[8]);
  res_t[2] =  F[0]*(F[2]*P[0] + F[5]*P[1] + F[8]*P[2]) + F[3]*(F[2]*P[3] + F[5]*P[4] + F[8]*P[5]) + F[6]*(F[2]*P[6] + F[5]*P[7] + F[8]*P[8]);
  // Second Column
  res_t[3] =  F[1]*(F[0]*P[0] + F[3]*P[1] + F[6]*P[2]) + F[4]*(F[0]*P[3] + F[3]*P[4] + F[6]*P[5]) + F[7]*(F[0]*P[6] + F[3]*P[7] + F[6]*P[8]);
  res_t[4] =  F[1]*(F[1]*P[0] + F[4]*P[1] + F[7]*P[2]) + F[4]*(F[1]*P[3] + F[4]*P[4] + F[7]*P[5]) + F[7]*(F[1]*P[6] + F[4]*P[7] + F[7]*P[8]);
  res_t[5] =  F[1]*(F[2]*P[0] + F[5]*P[1] + F[8]*P[2]) + F[4]*(F[2]*P[3] + F[5]*P[4] + F[8]*P[5]) + F[7]*(F[2]*P[6] + F[5]*P[7] + F[8]*P[8]);
  // Third Column
  res_t[6] =  F[2]*(F[0]*P[0] + F[3]*P[1] + F[6]*P[2]) + F[5]*(F[0]*P[3] + F[3]*P[4] + F[6]*P[5]) + F[8]*(F[0]*P[6] + F[3]*P[7] + F[6]*P[8]);
  res_t[7] =  F[2]*(F[1]*P[0] + F[4]*P[1] + F[7]*P[2]) + F[5]*(F[1]*P[3] + F[4]*P[4] + F[7]*P[5]) + F[8]*(F[1]*P[6] + F[4]*P[7] + F[7]*P[8]);
  res_t[8] =  F[2]*(F[2]*P[0] + F[5]*P[1] + F[8]*P[2]) + F[5]*(F[2]*P[3] + F[5]*P[4] + F[8]*P[5]) + F[8]*(F[2]*P[6] + F[5]*P[7] + F[8]*P[8]);

  res[0] = res_t[0];
  res[1] = res_t[1];
  res[2] = res_t[2];
  res[3] = res_t[3];
  res[4] = res_t[4];
  res[5] = res_t[5];
  res[6] = res_t[6];
  res[7] = res_t[7];
  res[8] = res_t[8];

}

void fromqtoR(float q[4], float R[9])
{
  // First column
  R[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  R[1] = 2*q[0]*q[3] + 2*q[1]*q[2];
  R[2] = 2*q[1]*q[3] - 2*q[0]*q[2];
  // Second column
  R[3] = 2*q[1]*q[2] - 2*q[0]*q[3];
  R[4] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  R[5] = 2*q[0]*q[1] + 2*q[2]*q[3];
  // Third column
  R[6] = 2*q[0]*q[2] + 2*q[1]*q[3];
  R[7] = 2*q[2]*q[3] - 2*q[0]*q[1];
  R[8] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

void leftQuaternion(float q[4], float qL[16])
{
  qL[0] =  q[0];
  qL[1] =  q[1];
  qL[2] =  q[2];
  qL[3] =  q[3];
  
  qL[4] = -q[1];
  qL[5] =  q[0];
  qL[6] =  q[3];
  qL[7] = -q[2];
  
  qL[8]  = -q[2];
  qL[9]  = -q[3];
  qL[10] =  q[0];
  qL[11] =  q[1];
  
  qL[12] = -q[3];
  qL[13] =  q[2];
  qL[14] = -q[1];
  qL[15] =  q[0];
  
}

void skew(float vec[3], float res[9])
{
  res[0] =  0.0;
  res[1] =  vec[2];
  res[2] = -vec[1];
  
  res[3] = -vec[2];
  res[4] =  0.0;
  res[5] =  vec[0];
  
  res[6] =  vec[1];
  res[7] = -vec[0];
  res[8] =  0.0;
}

void computeEulerAngles(float q[4], float euler[3])
{
    euler[0] = atan2(2.0f*(q[0]*q[1]+q[2]*q[3]),q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]);
    euler[1] = asin(2.0f*(q[1]*q[3]-q[0]*q[2]));
    euler[2] = atan2(2.0f*(q[1]*q[2]+q[0]*q[3]),q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]);
}


void matrixByVector3(float mat[9], float vec[3], float res[3])
{
  float res_t[3];
  res_t[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
  res_t[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
  res_t[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];
  
  for (int jj = 0; jj < 3; ++jj)
  {
    res[jj] = res_t[jj];
  }
}

void matrixByVector4(const float mat[16], float vec[4], float res[4])
{
    
  float res_t[4];
  res_t[0] = mat[0]*vec[0] + mat[4]*vec[1] + mat[8]*vec[2]  + mat[12]*vec[3];
  res_t[1] = mat[1]*vec[0] + mat[5]*vec[1] + mat[9]*vec[2]  + mat[13]*vec[3];
  res_t[2] = mat[2]*vec[0] + mat[6]*vec[1] + mat[10]*vec[2] + mat[14]*vec[3];
  res_t[3] = mat[3]*vec[0] + mat[7]*vec[1] + mat[11]*vec[2] + mat[15]*vec[3];

  for (int jj = 0; jj < 4; ++jj)
  {
    res[jj] = res_t[jj];  
  }
}

void matrixByMatrix(uint8_t dim, float matA[9], float matB[9], float res[9])
{
  // Multiplication of square matrices, res = matA*matB
  float res_t[dim*dim];
  for (int ii = 0; ii < dim; ++ii)
  {
    for (int jj = 0; jj < dim; ++jj)
    {
      for (int kk = 0; kk < dim; ++kk)
      {
        res_t[ii+dim*jj] += matA[ii+dim*kk]*matB[kk+dim*jj];
      }
    }
  }

  for (int ii = 0; ii < dim*dim; ++ii)
  {
    res[ii] = res_t[ii];
  }
}


void normalizeVector3(float vec[3], float res[3])
{
  float norm;
  norm = sqrtf(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

  res[0] = vec[0]/norm;
  res[1] = vec[1]/norm;
  res[2] = vec[2]/norm;
    
}

void normalizeVector4(float vec[4], float res[4])
{
  float norm;
  norm = sqrtf(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3]);
  
  res[0] = vec[0]/norm;
  res[1] = vec[1]/norm;
  res[2] = vec[2]/norm;
  res[3] = vec[3]/norm;
}


void addMatrices(uint8_t nrow, uint8_t ncol, float matA[9], float matB[9], float res[9])
{
  float res_t[nrow*ncol];
  for (int ii = 0; ii < nrow*ncol; ++ii)
  {
    res_t[ii] = matA[ii] + matB[ii];
  }

  for (int ii = 0; ii < nrow*ncol; ++ii)
  {
    res[ii] = res_t[ii];
  }
}

void subtractMatrices(uint8_t nrow, uint8_t ncol, float matA[9], float matB[9], float res[9])
{
  float res_t[nrow*ncol];
  for (int ii = 0; ii < nrow*ncol; ++ii)
  {
    res_t[ii] = matA[ii] - matB[ii];
  }
  for (int ii = 0; ii < nrow*ncol; ++ii)
  {
    res[ii] = res_t[ii];
  }
}

void inverseMatrix3(float mat[9], float res[9])
{
  float det = mat[0]*mat[4]*mat[8] - mat[0]*mat[5]*mat[7] - mat[1]*mat[3]*mat[8] + mat[1]*mat[5]*mat[6] + mat[2]*mat[3]*mat[7] - mat[2]*mat[4]*mat[6];

  // First column
  res[0] =  (mat[4]*mat[8] - mat[5]*mat[7])/det;
  res[1] = -(mat[1]*mat[8] - mat[2]*mat[7])/det;
  res[2] =  (mat[1]*mat[5] - mat[2]*mat[4])/det;
  // Second column
  res[3] = -(mat[3]*mat[8] - mat[5]*mat[6])/det;
  res[4] =  (mat[0]*mat[8] - mat[2]*mat[6])/det;
  res[5] = -(mat[0]*mat[5] - mat[2]*mat[3])/det;
  // Third column
  res[6] =  (mat[3]*mat[7] - mat[4]*mat[6])/det;
  res[7] = -(mat[0]*mat[7] - mat[1]*mat[6])/det;
  res[8] =  (mat[0]*mat[4] - mat[1]*mat[3])/det;
}

void transpose3(float mat[9], float res[9])
{
  res[0] = mat[0];
  res[1] = mat[3];
  res[2] = mat[6];
  res[3] = mat[1];
  res[4] = mat[4];
  res[5] = mat[7];
  res[6] = mat[2];
  res[7] = mat[5];
  res[8] = mat[8];
}

void printVector3(const char title[50], float vec[3])
{
  Serial.println(title);
  Serial.print("[ ");
  Serial.print(vec[0],12);
  Serial.print("  ");
  Serial.print(vec[1],12);
  Serial.print("  ");
  Serial.print(vec[2],12);  
  Serial.println(" ]'");
}

void printVector4(const char title[50], float vec[4])
{
  Serial.println(title);
  Serial.print("[ ");
  Serial.print(vec[0],12);
  Serial.print("  ");
  Serial.print(vec[1],12);
  Serial.print("  ");
  Serial.print(vec[2],12);
  Serial.print("  ");
  Serial.print(vec[3],12);    
  Serial.println(" ]'");
}

void printMatrix3(const char title[50], float mat[9])
{
  Serial.println(title);
  Serial.print("[[ ");
  Serial.print(mat[0],32);
  Serial.print("  ");
  Serial.print(mat[3],32);
  Serial.print("  ");
  Serial.print(mat[6],32);
  Serial.print("];[");
  Serial.print(mat[1],32);
  Serial.print("  ");
  Serial.print(mat[4],32);
  Serial.print("  ");
  Serial.print(mat[7],32);
  Serial.print("];[");
  Serial.print(mat[2],32);
  Serial.print("  ");
  Serial.print(mat[5],32);
  Serial.print("  ");
  Serial.print(mat[8],32);    
  Serial.println(" ]]");
}

void printMatrix4(const char title[50], float mat[16])
{
  Serial.println(title);
  Serial.print("[[ ");
  Serial.print(mat[0],12);
  Serial.print("  ");
  Serial.print(mat[4],12);
  Serial.print("  ");
  Serial.print(mat[8],12);
  Serial.print("  ");
  Serial.print(mat[12],12);
  Serial.print("];[");
  Serial.print(mat[1],12);
  Serial.print("  ");
  Serial.print(mat[5],12);
  Serial.print("  ");
  Serial.print(mat[9],12);
  Serial.print("  ");
  Serial.print(mat[13],12);
  Serial.print("];[");
  Serial.print(mat[2],12);
  Serial.print("  ");
  Serial.print(mat[6],12);
  Serial.print("  ");
  Serial.print(mat[10],12);
  Serial.print("  ");
  Serial.print(mat[14],12);
  Serial.print("];[");
  Serial.print(mat[3],12);
  Serial.print("  ");
  Serial.print(mat[7],12);
  Serial.print("  ");
  Serial.print(mat[11],12);
  Serial.print("  ");
  Serial.print(mat[15],12);      
  Serial.println(" ]]");
}


