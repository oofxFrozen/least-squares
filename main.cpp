// Nikita Sannikov

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

typedef vector<vector<double>> Dimensions;
typedef vector<double> Row;

class DimensionsException : public exception
{
public:
    DimensionsException() : msg("Error: the dimensional problem occurred") {}

    const char *what() const noexcept
    {
        return msg.c_str();
    }

private:
    string msg;
};

class IllegalOperationException : public exception
{
public:
    IllegalOperationException() : msg("Error: this operation cannot be applied to the matrix") {}

    const char *what() const noexcept
    {
        return msg.c_str();
    }

private:
    string msg;
};

class Matrix
{
private:
public:
    int n, m;
    Dimensions data;
    Matrix(int n, int m, Dimensions data)
    {
        this->n = n;
        this->m = m;
        this->data = data;
    }
    Matrix()
    {
    }
    Matrix operator+(const Matrix &other) const
    {
        if (this->n != other.n || this->m != other.m)
        {
            throw DimensionsException();
        }

        Matrix res(this->n, this->m, this->data);
        for (int i = 0; i < this->n; ++i)
            for (int j = 0; j < this->m; ++j)
            {
                res.data.at(i).at(j) += other.data.at(i).at(j);
            }
        return res;
    }
    Matrix operator-(const Matrix &other) const
    {
        if (this->n != other.n || this->m != other.m)
        {
            throw DimensionsException();
        }
        Matrix res(this->n, this->m, this->data);
        for (int i = 0; i < this->n; ++i)
            for (int j = 0; j < this->m; ++j)
            {
                res.data.at(i).at(j) -= other.data.at(i).at(j);
            }
        return res;
    }
    void operator=(const Matrix &other)
    {
        setN(other.n);
        setM(other.m);
        setDimensions(other.data);
    }
    Row operator[](int i) const
    {
        return data[i];
    }
    friend istream &operator>>(istream &in, Matrix &mat)
    {
        int n, m;
        double temp;

        in >> n >> m;
        Dimensions data = Dimensions(n, Row(m, 0));

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                in >> temp;
                data.at(i).at(j) = temp;
            }
        }

        mat = Matrix(n, m, data);

        return in;
    }
    friend void operator<<(ostream &out, const Matrix &m)
    {
        cout.precision(4);
        cout << fixed;
        for (Row row : m.data)
        {
            for (double a : row)
                out << a << " ";
            out << endl;
        }
    }
    Matrix operator*(const Matrix &other) const
    {
        if (this->m != other.n)
        {
            throw DimensionsException();
        }

        Matrix res(this->n, other.m, Dimensions());
        for (int i = 0; i < this->n; ++i)
        {
            res.data.push_back(Row());
            for (int j = 0; j < other.m; ++j)
            {
                res.data.at(i).push_back(0);
                double sum = 0;
                for (int l = 0; l < this->m; ++l)
                    sum += this->data.at(i).at(l) * other.data.at(l).at(j);
                res.data.at(i).at(j) = sum;
            }
        }
        return res;
    }
    Matrix transposed() const
    {
        Dimensions dim = Dimensions();
        for (int i = 0; i < m; i++)
        {
            Row row = Row();
            for (int j = 0; j < n; j++)
            {
                row.push_back(data.at(j).at(i));
            }
            dim.push_back(row);
        }

        return Matrix(m, n, dim);
    }
    void swapRows(int n, int m)
    {
        Row b = data[n];
        data[n] = data[m];
        data[m] = b;
    }
    void setN(int newN)
    {
        this->n = newN;
    }
    void setM(int newM)
    {
        this->m = newM;
    }
    void setDimensions(Dimensions newDimensions)
    {
        this->data = newDimensions;
    }
};

class DataMatrix : public Matrix
{
private:
public:
    DataMatrix(int n, Dimensions data) : Matrix(n, 2, data)
    {
    }
    DataMatrix() : Matrix() {}
    friend istream &operator>>(istream &in, DataMatrix &mat)
    {
        int n;
        double temp;

        in >> n;
        Dimensions data = Dimensions(n, Row(2, 0));

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                in >> temp;
                data.at(i).at(j) = temp;
            }
        }
        mat = DataMatrix(n, data);

        return in;
    }
};

class SquareMatrix : public Matrix
{
public:
    SquareMatrix(int n, Dimensions data) : Matrix(n, n, data)
    {
    }
    SquareMatrix() : Matrix() {}
    SquareMatrix diagonalMatrix()
    {
        Dimensions data = Dimensions(n, Row(n, 0));
        for (int i = 0; i < n; ++i)
        {
            data.at(i).at(i) = this->data[i][i];
        }
        return SquareMatrix(n, data);
    }
    friend istream &operator>>(istream &in, SquareMatrix &mat)
    {
        int n;
        int temp;

        in >> n;
        Dimensions data = Dimensions(n, Row(n, 0));

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                in >> temp;
                data.at(i).at(j) = temp;
            }
        }
        mat = SquareMatrix(n, data);

        return in;
    }
    SquareMatrix operator*(const SquareMatrix &other) const
    {
        if (this->m != other.n)
        {
            throw DimensionsException();
        }

        SquareMatrix res(this->n, Dimensions());
        for (int i = 0; i < this->n; ++i)
        {
            res.data.push_back(Row());
            for (int j = 0; j < other.m; ++j)
            {
                double sum = 0;
                for (int l = 0; l < this->m; ++l)
                    sum += this->data.at(i).at(l) * other.data.at(l).at(j);
                res.data.at(i).push_back(sum);
            }
        }
        return res;
    }
    Matrix operator*(const Matrix &other) const
    {
        if (this->m != other.n)
        {
            throw DimensionsException();
        }
        Matrix res(this->n, other.m, Dimensions());
        for (int i = 0; i < this->n; ++i)
        {
            res.data.push_back(Row());
            for (int j = 0; j < other.m; ++j)
            {
                double sum = 0;
                for (int l = 0; l < this->m; ++l)
                    sum += this->data.at(i).at(l) * other.data.at(l).at(j);
                res.data.at(i).push_back(sum);
            }
        }
        return res;
    }
    SquareMatrix operator-(const Matrix &other) const
    {
        if (this->n != other.n || this->m != other.m)
        {
            throw DimensionsException();
        }

        SquareMatrix res(this->n, this->data);
        for (int i = 0; i < this->n; ++i)
            for (int j = 0; j < this->m; ++j)
            {
                res.data.at(i).at(j) -= other.data.at(i).at(j);
            }
        return res;
    }
    SquareMatrix lower_triangular()
    {
        Dimensions data = this->data;
        for (int i = 0; i < this->n - 1; ++i)
        {
            for (int j = i + 1; j < this->n; ++j)
            {
                data.at(i).at(j) = 0;
            }
        }
        return SquareMatrix(this->n, data);
    }
    SquareMatrix upper_triangular()
    {
        Dimensions data = this->data;
        for (int i = 0; i < this->n; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                data.at(i).at(j) = 0;
            }
        }
        return SquareMatrix(this->n, data);
    }
};

class IdentityMatrix : public SquareMatrix
{
public:
    IdentityMatrix(int n)
    {
        Row r;
        Dimensions data = Dimensions();
        for (int i = 0; i < n; i++)
        {
            r = Row();
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    r.push_back(1);
                }
                else
                {
                    r.push_back(0);
                }
            }
            data.push_back(r);
        }

        this->SquareMatrix::~SquareMatrix();
        new (this) SquareMatrix(n, data);
    }

    friend void operator>>(istream &in, IdentityMatrix &mat)
    {
        throw IllegalOperationException();
    }
};

class EliminationMatrix : public SquareMatrix
{
public:
    EliminationMatrix(SquareMatrix &initial, int n, int m)
    {
        Row r;
        Dimensions data = Dimensions();
        for (int i = 0; i < initial.n; i++)
        {
            r = Row();
            for (int j = 0; j < initial.n; j++)
            {
                if (i == n && j == m)
                {
                    r.push_back(-initial[n][m] / initial[m][m]);
                    continue;
                }
                if (i == j)
                {
                    r.push_back(1);
                }
                else
                {
                    r.push_back(0);
                }
            }
            data.push_back(r);
        }

        this->SquareMatrix::~SquareMatrix();
        new (this) SquareMatrix(initial.n, data);
    }
    friend void operator>>(istream &in, EliminationMatrix &mat)
    {
        throw IllegalOperationException();
    }
};

class PermutationMatrix : public SquareMatrix
{
public:
    PermutationMatrix(SquareMatrix &initial, int n, int m)
    {
        IdentityMatrix i = IdentityMatrix(initial.n);
        i.swapRows(n, m);

        this->SquareMatrix::~SquareMatrix();
        new (this) SquareMatrix(i.n, i.data);
    }
    friend void operator>>(istream &in, PermutationMatrix &mat)
    {
        throw IllegalOperationException();
    }
};

SquareMatrix inverse_matrix(SquareMatrix matrix)
{
    SquareMatrix mat = matrix;
    SquareMatrix identity = SquareMatrix(matrix.n, IdentityMatrix(matrix.n).data);
    int step = 1;

    for (int i = 0; i < mat.n; ++i)
    {
        double max = i;
        for (int j = i + 1; j < mat.n; ++j)
        {
            if (abs(mat[j][i]) > abs(mat[max][i]))
            {
                max = j;
            }
        }
        if (max != i)
        {
            mat = PermutationMatrix(mat, max, i) * mat;
            identity = PermutationMatrix(mat, max, i) * identity;
        }

        for (int j = i + 1; j < mat.n; ++j)
        {
            if (mat[j][i] != 0)
            {
                EliminationMatrix elim = EliminationMatrix(mat, j, i);
                mat = elim * mat;
                identity = elim * identity;
            }
        }
    }

    for (int i = mat.n - 1; i >= 0; --i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            if (mat[j][i] != 0)
            {
                EliminationMatrix elim = EliminationMatrix(mat, j, i);
                mat = elim * mat;
                identity = elim * identity;
            }
        }
    }

    for (int i = 0; i < mat.n; ++i)
    {
        double coeff = mat[i][i];
        mat.data.at(i).at(i) = 1;

        for (int j = 0; j < mat.n; ++j)
        {
            identity.data.at(i).at(j) = identity[i][j] / coeff;
        }
    }

    return identity;
}

Matrix getA(Matrix data, int degree)
{
    Matrix A = Matrix(data.n, degree + 1, Dimensions());
    for (int i = 0; i < data.n; ++i)
    {
        Row r = Row();
        double power = 1;
        for (int j = 0; j < degree + 1; ++j)
        {
            r.push_back(power);
            power *= data.data[i][0];
        }
        A.data.push_back(r);
    }
    return A;
}

Matrix getB(Matrix data)
{
    Matrix b = Matrix(data.n, 1, Dimensions());
    for (int i = 0; i < data.n; ++i)
    {
        Row r = Row();
        r.push_back(data.data.at(i).at(1));
        b.data.push_back(r);
    }
    return b;
}

double evaluate_curve(double x, const vector<double> &coeffs)
{
    double y = 0.0;
    double pow = 1;
    for (int i = 0; i < coeffs.size(); i++)
    {
        y += coeffs[i] * pow;
        pow *= x;
    }
    return y;
}

vector<double> least_squares(DataMatrix data, int degree)
{
    Matrix A = getA(data, degree);
    Matrix b = getB(data);
    Matrix ATA = A.transposed() * A;
    SquareMatrix ATA_sqr = SquareMatrix(ATA.n, ATA.data);
    SquareMatrix ATA_1 = inverse_matrix(ATA_sqr);

    cout
        << "A:\n"
        << A;
    cout
        << "A_T*A:\n"
        << ATA;

    cout
        << "(A_T*A)^-1:\n"
        << ATA_1;

    cout
        << "A_T*b:\n"
        << A.transposed() * b;
    Matrix result = ATA_1 * A.transposed() * b;
    cout
        << "x~:\n"
        << result;

    vector<double> coefficients = vector<double>(result.n);
    for (int i = 0; i < result.n; ++i)
    {
        coefficients.at(i) = result[i][0];
    }

    return coefficients;
}

main()
{

    DataMatrix data;
    int degree;
    cin >> data >> degree;

    vector<double> coefficients = least_squares(data, degree);

    // Generate points on the curve defined by the coefficients
    std::ofstream outfile("results.dat");
    for (double xi = -20; xi <= 20; xi += 0.1)
    {
        double yi = evaluate_curve(xi, coefficients);
        outfile << xi << " " << yi << std::endl;
    }
    outfile.close();

    return 0;
}
