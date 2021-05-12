#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

//矩阵类的定义
class matrix {
public:
    int  m_Row; // 矩阵的行数
    int  m_Col; // 矩阵的列数
    double* m_lpBuf; // 动态分配用来存放数组的空间
public:
    matrix(int m, int n);
    matrix(const matrix& A);
    ~matrix();
    matrix(int rows, int cols, double value);

    //矩阵乘法
    matrix multiply(matrix& a, matrix& b) {

        if (a.m_Col != b.m_Row) {
            cout << "无法相乘" << endl;
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

    //矩阵求逆矩阵
    matrix inv(matrix a) {
        int rows = a.m_Row;
        matrix C(rows, rows);
        if (a.m_Col != a.m_Row) {
            cout << "必须是方阵" << endl;
            return C;
        }

        // 构造一个与A行相同的单位阵B = |A : E|
        matrix B(rows, rows * 2, 0);
        for (int i = 1; i <= rows; ++i)
            for (int j = 1; j <= rows; ++j)
                B(i, j) = a(i, j);
        for (int i = 1; i <= rows; ++i)
            B(i, i + rows) = 1;

        double tempa = 0;
        // 对矩阵A进行B.row次消元运算，每次保证第K列只有对角线上非零 
        for (int k = 1; k <= B.m_Row; k++)
        {

            double max = fabs(B(k, k)); // 主元初始默认为右下方矩阵首个元素 
            int ind = k; // 主元行号默认为右下方矩阵首行 
            // 结果第ind行为列主元行 
            for (int n = k + 1; n <= B.m_Row; ++n)
            {
                if (fabs(B(n, k)) > max)
                {
                    // 遇到绝对值更大的元素 
                    max = fabs(B(n, k)); // 更新主元值 
                    ind = n; // 更新主元行号
                }
            }

            if (ind != k)
            { // 主元行不是右下方矩阵首行 
                for (int m = k; m <= rows * 2; ++m)
                { // 将主元行与右下方矩阵首行互换 
                    tempa = B(k, m);
                    B(k, m) = B(ind, m);
                    B(ind, m) = tempa;
                }
            }

            // 第k次消元操作，以第k行作为主元行，将其上下各行的第k列元素化为零 
            // 同时以同样的参数对B施以同样的操作，此时可以将B看作A矩阵的一部分 
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

    double& operator ()(int i, int j)//重载（）
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
matrix::matrix(int m, int n)//初始化矩阵
{
    m_Row = m;
    m_Col = n;
    this->m_lpBuf = new double[m * n];
    memset(m_lpBuf, 0, m * n * sizeof(double));
}
matrix& matrix::operator = (const matrix& A) //重载赋值运算符
{
    if (this == &A)
        return *this;
    this->resize(A.m_Row, A.m_Col);

    for (int i = 0; i < (m_Col * m_Row); i++)
        m_lpBuf[i] = A.m_lpBuf[i];
    return *this;
}
    
matrix::matrix(const matrix& A) //拷贝构造函数
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
matrix::matrix(int rows, int cols, double value)//初始化为特殊值
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
    infile.open(pathname.data());   //将文件流对象与文件连接起来
    assert(infile.is_open());   //若失败,则输出错误消息,并终止程序运行
    vector<double> suanz;
    string s;
    while (getline(infile, s)) {
        if (s[0] == '#')
            continue;
        istringstream is(s); //将读出的一行转成数据流进行操作
        double d;
        while (!is.eof()) {
            is >> d;
            suanz.push_back(d);
        }
        res.push_back(suanz);
        suanz.clear();
        s.clear();
    }
    infile.close();             //关闭文件输入流
}

int main() {
    vector <vector<double>> readIn;
    vector <double> innerElement;   // 内方位元素
    vector <double> m;  // 摄影比例尺分母
    vector <double> controlPointsNumber;    // 控制点对数
    vector <vector<double>> points; // 控制点对 x, y, X, Y, Z

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


    //确定未知数的初始值
    double Xs = (l_x[0] + l_x[1] + l_x[2] + l_x[3]) / 4;
    double Ys = (l_y[0] + l_y[1] + l_y[2] + l_y[3]) / 4;
    double Zs = f * m[0];
    double φ0, ω0, κ0 = 0;

    //计算旋转矩阵
    matrix R(3,3,0);
}

void Rotation_matrix(matrix R, double φ0, double ω0, double κ0) {
    R(1,1) = cos(φ0) * cos(κ0) - sin(φ0) * sin(ω0) * sin(κ0);
    R(1,2)= -cos(φ0) * sin(κ0) - sin(φ0) * sin(ω0) * cos(κ0);
    R(1,3)=  - sin(φ0) * cos(ω0) ;
    R(2,1)=  cos(ω0) * sin(κ0);
    R(2,2)= cos(ω0) * cos(κ0);
    R(2,3)= -sin(ω0) ;
    R(3,1)= sin(φ0) * cos(κ0) + cos(φ0) * sin(ω0) * sin(κ0);
    R(3,2) = -sin(φ0) * sin(κ0) + cos(φ0) * sin(ω0) * cos(κ0);
    R(3,3)= cos(φ0) * cos(ω0) ;
}

double ApproximateValueX(double x0,double f, matrix R,double Xs,double Ys,double Zs,double X,double Y,double Z) {
    return x0 - f * (R(1,1) * (X - Xs) + R(2,1) * (Y - Ys) + R(3,1) * (Z - Zs)) / (R(1,3) * (X - Xs) + R(2,3) * (Y - Ys) + R(3,3) * (Z - Zs));
}

double ApproximateValueY(double y0, double f, matrix R, double Xs, double Ys, double Zs, double X, double Y, double Z) {
    return y0 - f * (R(1, 2) * (X - Xs) + R(2, 2) * (Y - Ys) + R(3, 2) * (Z - Zs)) / (R(1, 3) * (X - Xs) + R(2, 3) * (Y - Ys) + R(3, 3) * (Z - Zs));
}




