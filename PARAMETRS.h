#pragma once

//----------------��������������� ���������---------------//
#define PI 3.1415926535897932384626433832795029

//----------------------�����-�����-----------------------//
#define M 10000   // ����� �������� ��� �����-�����


//-------------------���������� �������-------------------//
#define x1 0.
#define x2 10.
#define L (x2-x1)

#define t1 0.
#define t2 10.

//------------------------�����---------------------------//
#define h 0.02 // ��� �� ����������
#define t 0.01 // ��� �� �������
#define r (t/h)

#define N ((x2 - x1) / h)
#define T ((t2 - t1) / t)

#define SIZE_X (int)(N)
#define SIZE_T (int)(T)


//-----------������������ ��������� ������-�������--------//
#define m 2.
#define g 10. // 1, 10
#define g_0 1.

//--------------------������������ ������-----------------//
#define TMP 10.  // ����������� (�����, ����� T > m)


//---------------------������������ �����-----------------//
#define al 10.      // 
#define x_q 5.    // ����� ����������
#define sigma ((10.*h)/(2*sqrt(2*log(2)))) //  = 0.085


