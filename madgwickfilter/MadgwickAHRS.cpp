
#include "MadgwickAHRS.h"
#include <math.h>

#define sampleFreqDef   512.0f          // sample frequency in Hz
#define betaDef         0.1f    
#define gyroMeasError 3.14159265358979 * (5.0f / 180.0f)
#define gyroMeasDrift 3.14159265358979 * (0.2f / 180.0f)


Madgwick::Madgwick() {
	beta = sqrt(3.0f / 4.0f) * gyroMeasError;
	q0 = 1.0f;
	q1 = 0.0f;
	q2 = 0.0f;
	q3 = 0.0f;
	float deltat = 1.0f / sampleFreqDef;
	anglesComputed = 0;
	
}

void Madgwick::updateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
	float beta = sqrt(3.0f / 4.0f) * gyroMeasError;
	float zeta = sqrt(3.0f / 4.0f) * gyroMeasDrift;
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;
	float b_x = 1, b_z = 0;
	float g_bx = 0, g_by = 0, g_bz = 0;
	float m_x,m_y,m_z;
	float g_err_x,g_err_y,g_err_z;
	float deltat = 1.0f / sampleFreqDef;
	
	//-------definations
	float f_1, f_2, f_3, f_4, f_5, f_6; // objective function elements
	float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33, J_41, J_42, J_43, J_44, J_51, J_52, J_53, J_54, J_61, J_62, J_63, J_64; //
    float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4; // estimated direction of the gyroscope error
float w_err_x, w_err_y, w_err_z; // estimated direction of the gyroscope error (angular)
float h_x, h_y, h_z; // computed flux in the earth frame
// axulirary variables to avoid reapeated calcualtions
float halfq0 = 0.5f * q0;
float halfq1 = 0.5f * q1;
float halfq2 = 0.5f * q2;
float halfq3 = 0.5f * q3;
float twoq0 = 2.0f * q0;
float twoq1 = 2.0f * q1;
float twoq2 = 2.0f * q2;
float twoq3= 2.0f * q3;
float twob_x = 2.0f * b_x;
float twob_z = 2.0f * b_z;
float twob_xq0 = 2.0f * b_x * q0;
float twob_xq1 = 2.0f * b_x * q1;
float twob_xq2 = 2.0f * b_x * q2;
float twob_xq3 = 2.0f * b_x * q3;
float twob_zq0 = 2.0f * b_z * q0;
float twob_zq1= 2.0f * b_z * q1;
float twob_zq2= 2.0f * b_z * q2;
float twob_zq3 = 2.0f * b_z * q3;
float q0q1=q0*q1;;
float q0q2 = q0 *q2;
float q0q3=q0*q3;
float q1q2=q1*q2;
float q1q3= q1*q3;
float q2q3=q2*q3;
float twom_x = 2.0f * m_x;
float twom_y = 2.0f * m_y;
float twom_z = 2.0f * m_z;
	//--------jacobians ahead
f_1 = twoq1 * q3 - twoq0 * q2 - ax;
f_2 = twoq0 * q1 + twoq2 * q3 - ay;
f_3 = 1.0f - twoq1 * q1 - twoq2 * q2 - az;
f_4 = twob_x * (0.5f - q2 * q2 - q3 * q3) + twob_z * (q1q3 - q0q2) - m_x;
f_5 = twob_x * (q1 * q2 - q0 * q3) + twob_z * (q0 * q1 + q2 * q3) - m_y;
f_6 = twob_x * (q0q2 + q1q3) + twob_z * (0.5f - q1 * q1 - q2 * q2) - m_z;
J_11or24 = twoq2; // J_11 negated in matrix multiplication
J_12or23 = 2.0f * q3;
J_13or22 = twoq0; // J_12 negated in matrix multiplication
J_14or21 = twoq1;
J_32 = 2.0f * J_14or21; // negated in matrix multiplication
J_33 = 2.0f * J_11or24; // negated in matrix multiplication
J_41 = twob_zq2; // negated in matrix multiplication
J_42 = twob_zq3;
J_43 = 2.0f * twob_xq2 + twob_zq0; // negated in matrix multiplication
J_44 = 2.0f * twob_xq3 - twob_zq1; // negated in matrix multiplication
J_51 = twob_xq3 - twob_zq1; // negated in matrix multiplication
J_52 = twob_xq2 + twob_zq0;
J_53 = twob_xq1+ twob_zq3;
J_54 = twob_xq1 - twob_zq2; // negated in matrix multiplication
J_61 = twob_xq2;
J_62 = twob_xq3 - 2.0f * twob_zq1;
J_63 = twob_xq0- 2.0f * twob_zq2;
J_64 = twob_xq1;
//gradients	
SEqHatDot_1 = J_14or21 * f_2 - J_11or24 * f_1 - J_41 * f_4 - J_51 * f_5 + J_61 * f_6;
SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2 - J_32 * f_3 + J_42 * f_4 + J_52 * f_5 + J_62 * f_6;
SEqHatDot_3 = J_12or23 * f_2 - J_33 * f_3 - J_13or22 * f_1 - J_43 * f_4 + J_53 * f_5 + J_63 * f_6;
SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2 - J_44 * f_4 - J_54 * f_5 + J_64 * f_6;
float norm = sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
SEqHatDot_1 = SEqHatDot_1 / norm;
SEqHatDot_2 = SEqHatDot_2 / norm;
SEqHatDot_3 = SEqHatDot_3 / norm;
SEqHatDot_4 = SEqHatDot_4 / norm;	

	
	
	//--------
	    g_err_x = twoq0 * SEqHatDot_2 - twoq1 * SEqHatDot_1 - twoq2 * SEqHatDot_4 + twoq3 * SEqHatDot_3;
        g_err_y = twoq0 * SEqHatDot_3 + twoq1 * SEqHatDot_4 - twoq2 * SEqHatDot_1 - twoq3* SEqHatDot_2;
        g_err_z = twoq0 * SEqHatDot_4 - twoq1 * SEqHatDot_3 + twoq2 * SEqHatDot_2 - twoq3 * SEqHatDot_1;
	gx *= 0.0174533f;
	gy *= 0.0174533f;
	gz *= 0.0174533f;
	g_bx += g_err_x * deltat * zeta;
        g_by += g_err_y * deltat * zeta;
        g_bz += g_err_z * deltat * zeta;
        gx -= g_bx;
        gy -= g_by;
        gz -= g_bz; 


	qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;   

		
		_2q0 = 2.0f * q0;
		_2q1 = 2.0f * q1;
		_2q2 = 2.0f * q2;
		_2q3 = 2.0f * q3;
		_4q0 = 4.0f * q0;
		_4q1 = 4.0f * q1;
		_4q2 = 4.0f * q2;
		_8q1 = 8.0f * q1;
		_8q2 = 8.0f * q2;
		q0q0 = q0 * q0;
		q1q1 = q1 * q1;
		q2q2 = q2 * q2;
		q3q3 = q3 * q3;

		
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
		float Norm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); 
		s0 *= Norm;
		s1 *= Norm;
		s2 *= Norm;
		s3 *= Norm;
	}

	
        q0 += qDot1 * deltat;
	q1 += qDot2 * deltat;
	q2 += qDot3 * deltat;
	q3 += qDot4 * deltat;

	
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
	anglesComputed = 0;
}

float Madgwick::invSqrt(float x) {
	float halfx = 0.5f * x;
	union { float f; long l; } i;
	i.f = x;
	i.l = 0x5f3759df - (i.l >> 1);
	float y = i.f;
	y = y * (1.5f - (halfx * y * y));
	y = y * (1.5f - (halfx * y * y));
	return y;
}

void Madgwick::computeAngles()
{
	roll = atan2f(q0*q1 + q2*q3, 0.5f - q1*q1 - q2*q2);
	pitch = asinf(-2.0f * (q1*q3 - q0*q2));
	yaw = atan2f(q1*q2 + q0*q3, 0.5f - q2*q2 - q3*q3);
	anglesComputed = 1;
}

