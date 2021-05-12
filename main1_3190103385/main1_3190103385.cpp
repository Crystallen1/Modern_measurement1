#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

//������Ķ���
class matrix {
public:
    int  m_Row; // ���������
    int  m_Col; // ���������
    double* m_lpBuf; // ��̬���������������Ŀռ�
public:
    matrix(int m, int n);
    matrix(const matrix& A);
    ~matrix();
    matrix(int rows, int cols, double value);

    //����˷�
    matrix multiply(matrix& a, matrix& b) {

        if (a.m_Col != b.m_Row) {
            cout << "�޷����" << endl;
            return;
        }

        matrix c(a.m_Row, b.m_Col,0);

        for (int m = 0; m < a.m_Row; m++) {
            for (int s = 0; s < b.m_Col; s++) {
                for (int n = 0; n < a.m_Col; n++) {
                    c(m,s)+= a(m,n) * b(n,s);
                }
            }   
        }
        return c;
    }

    //�����������
    matrix inv(matrix a) {
        int rows = a.m_Row;
        matrix C(rows, rows);
        if (a.m_Col != a.m_Row) {
            cout << "�����Ƿ���" << endl;
            return C;
        }

        // ����һ����A����ͬ�ĵ�λ��B = |A : E|
        matrix B(rows, rows * 2, 0);
        for (int i = 1; i <= rows; ++i)
            for (int j = 1; j <= rows; ++j)
                B(i, j) = a(i, j);
        for (int i = 1; i <= rows; ++i)
            B(i, i + rows) = 1;

        double tempa = 0;
        // �Ծ���A����B.row����Ԫ���㣬ÿ�α�֤��K��ֻ�жԽ����Ϸ��� 
        for (int k = 1; k <= B.m_Row; k++)
        {

            double max = fabs(B(k, k)); // ��Ԫ��ʼĬ��Ϊ���·������׸�Ԫ�� 
            int ind = k; // ��Ԫ�к�Ĭ��Ϊ���·��������� 
            // �����ind��Ϊ����Ԫ�� 
            for (int n = k + 1; n <= B.m_Row; ++n)
            {
                if (fabs(B(n, k)) > max)
                {
                    // ��������ֵ�����Ԫ�� 
                    max = fabs(B(n, k)); // ������Ԫֵ 
                    ind = n; // ������Ԫ�к�
                }
            }

            if (ind != k)
            { // ��Ԫ�в������·��������� 
                for (int m = k; m <= rows * 2; ++m)
                { // ����Ԫ�������·��������л��� 
                    tempa = B(k, m);
                    B(k, m) = B(ind, m);
                    B(ind, m) = tempa;
                }
            }

            // ��k����Ԫ�������Ե�k����Ϊ��Ԫ�У��������¸��еĵ�k��Ԫ�ػ�Ϊ�� 
            // ͬʱ��ͬ���Ĳ�����Bʩ��ͬ���Ĳ�������ʱ���Խ�B����A�����һ���� 
            tempa = 1.0 / B(k, k);
            for (int i = k; i <= rows * 2; ++i)
            {
                B(k, i) *= tempa;
            }

            for (int i = 1; i <= rows; i++)
            {
                if (i != k)
                {
                    tempa = -B(i, k);
                    for (int j = k; j <= B.m_Col; ++j)
                        B(i, j) += tempa * B(k, j);
                }
            }

        }

        for (int i = 1; i <= rows; ++i)
            for (int j = 1; j <= rows; ++j)
                C(i, j) = B(i, j + rows);

        return C;
    }

    double& operator ()(int i, int j)//���أ���
    {
        return *(m_lpBuf + (i - 1) * m_Col + (j - 1));
    }

    void resize(int row, int col) {
        m_Col = col;
        m_Row = row;
    }

    ~matrix() {
        m_Row = 0;
        m_Col = 0;
        if (m_lpBuf != NULL)
            delete[] this->m_lpBuf;
    };
};
matrix::matrix(int m, int n)//��ʼ������
{
    m_Row = m;
    m_Col = n;
    this->m_lpBuf = new double[m * n];
    memset(m_lpBuf, 0, m * n * sizeof(double));
}
matrix& matrix::operator = (const matrix& A) //���ظ�ֵ�����
{
    if (this == &A)
        return *this;
    this->resize(A.m_Row, A.m_Col);

    for (int i = 0; i < (m_Col * m_Row); i++)
        m_lpBuf[i] = A.m_lpBuf[i];
    return *this;
}
    
matrix::matrix(const matrix& A) //�������캯��
{
    this->m_Row = A.m_Row;
    this->m_Col = A.m_Col;
    if (m_Row * m_Col != 0)
    {
        this->m_lpBuf = new double[m_Row * m_Col];

        for (int i = 0; i < (m_Row * m_Col); ++i)
            this->m_lpBuf[i] = A.m_lpBuf[i];
    }
    else
        m_lpBuf = NULL;
}
matrix::matrix(int rows, int cols, double value)//��ʼ��Ϊ����ֵ
{
    m_Row = rows;
    m_Col = cols;
    for (int i = 0; i < m_Row; i++) {
        for (int j = 0; j < m_Col; j++) {
            m_lpBuf[i*m_Col+j] = value;
        }
    }
}

void txt_to_vector(vector<vector<double>>& res, string pathname)
{
    ifstream infile;
    infile.open(pathname.data());   //���ļ����������ļ���������
    assert(infile.is_open());   //��ʧ��,�����������Ϣ,����ֹ��������
    vector<double> suanz;
    string s;
    while (getline(infile, s)) {
        if (s[0] == '#')
            continue;
        istringstream is(s); //��������һ��ת�����������в���
        double d;
        while (!is.eof()) {
            is >> d;
            suanz.push_back(d);
        }
        res.push_back(suanz);
        suanz.clear();
        s.clear();
    }
    infile.close();             //�ر��ļ�������
}

int main() {
    vector <vector<double>> readIn;
    vector <double> innerElement;   // �ڷ�λԪ��
    vector <double> m;  // ��Ӱ�����߷�ĸ
    vector <double> controlPointsNumber;    // ���Ƶ����
    vector <vector<double>> points; // ���Ƶ�� x, y, X, Y, Z

    double x0, y0, f;

    txt_to_vector(readIn, "data1.txt");

    innerElement = readIn[0];
    m = readIn[1];
    controlPointsNumber = readIn[2];
    const int  controlPointsNumber1 = controlPointsNumber[0];
    double p_x[4], p_y[4], l_x[4], l_y[4], l_z[4];

    for (int i = 3; i < readIn.size(); i++) {
        points.push_back(readIn[i]);
    }

    x0=innerElement[0] ;
    y0=innerElement[1] ;
    f= innerElement[2] ;
    for (int i = 0; i <= 3; i++) {
        p_x[i] = points[i][0];
        p_y[i] = points[i][1];
        l_x[i] = points[i][2];
        l_y[i] = points[i][3];
        l_z[i] = points[i][4];
    }


    //ȷ��δ֪���ĳ�ʼֵ
    double Xs = (l_x[0] + l_x[1] + l_x[2] + l_x[3]) / 4;
    double Ys = (l_y[0] + l_y[1] + l_y[2] + l_y[3]) / 4;
    double Zs = f * m[0];
    double ��0, ��0, ��0 = 0;

    //������ת����
    matrix R(3,3,0);
}

void Rotation_matrix(matrix R, double ��0, double ��0, double ��0) {
    R(1,1) = cos(��0) * cos(��0) - sin(��0) * sin(��0) * sin(��0);
    R(1,2)= -cos(��0) * sin(��0) - sin(��0) * sin(��0) * cos(��0);
    R(1,3)=  - sin(��0) * cos(��0) ;
    R(2,1)=  cos(��0) * sin(��0);
    R(2,2)= cos(��0) * cos(��0);
    R(2,3)= -sin(��0) ;
    R(3,1)= sin(��0) * cos(��0) + cos(��0) * sin(��0) * sin(��0);
    R(3,2) = -sin(��0) * sin(��0) + cos(��0) * sin(��0) * cos(��0);
    R(3,3)= cos(��0) * cos(��0) ;
}

double ApproximateValueX(double x0,double f, matrix R,double Xs,double Ys,double Zs,double X,double Y,double Z) {
    return x0 - f * (R(1,1) * (X - Xs) + R(2,1) * (Y - Ys) + R(3,1) * (Z - Zs)) / (R(1,3) * (X - Xs) + R(2,3) * (Y - Ys) + R(3,3) * (Z - Zs));
}

double ApproximateValueY(double y0, double f, matrix R, double Xs, double Ys, double Zs, double X, double Y, double Z) {
    return y0 - f * (R(1, 2) * (X - Xs) + R(2, 2) * (Y - Ys) + R(3, 2) * (Z - Zs)) / (R(1, 3) * (X - Xs) + R(2, 3) * (Y - Ys) + R(3, 3) * (Z - Zs));
}




