#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

vector<vector<double>> extraiColunas(const vector<vector<double>> &A, const vector<int>& B) {
    vector<vector<double>> resultado(A.size(), vector<double>(B.size()));
    for (int i = 0; i < A.size(); ++i)
        for (int j = 0; j < B.size(); ++j)
            resultado[i][j] = A[i][B[j]];
    return resultado;
}

/**
 * Verifica se um vetor contém um valor específico.
 * @param v     [IN]  O vetor a ser verificado.
 * @param valor [IN]  O valor a ser procurado.
 * @retval bool  Retorna true se o valor estiver presente, false caso contrário.
 */
bool contemValor(vector<int> v, int valor)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == valor)
        {
            return true;
        }
    }
    return false;
}

/**
 * Verifica se um vetor contém um conjunto específico.
 * @param v     [IN]  O vetor a ser verificado.
 * @param conjunto [IN]  O conjunto a ser procurado.
 * @retval bool  Retorna true se o conjunto estiver presente, false caso contrário.
 */
bool contemConjunto(vector<vector<int>> v, vector<int> conjunto)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] == conjunto)
        {
            return true;
        }
    }
    return false;
}

/**
 * Encontra o menor valor em um vetor.
 * @param v     [IN]  O vetor a ser verificado.
 * @retval T  Retorna o menor valor encontrado.
 */
template <typename T>
T menorValor(vector<T> v)
{
    T menor = v[0];
    for (int i = 1; i < v.size(); i++)
    {
        if (v[i] < menor)
        {
            menor = v[i];
        }
    }
    return menor;
}

/**
 * Verifica se todos os elementos de um vetor são menores ou iguais a zero.
 * @param v     [IN]  O vetor a ser verificado.
 * @retval bool  Retorna true se todos os elementos forem menores ou iguais a zero, false caso contrário.
 */
template <typename T>
bool vetorMenorZero(vector<T> v)
{
    for (int i = 0; i < v.size(); i++)
    {
        if (v[i] > 0)
        {
            return false;
        }
    }
    return true;
}

template <typename T>
void imprimeMatriz(const vector<vector<T>> &M)
{
    for (const auto &row : M)
    {
        for (const auto &elem : row)
        {
            cout << '\t' << elem;
        }
        cout << endl;
    }
}

template <typename T>
void impremeVetor(const vector<T> &v)
{
    for (const auto &elem : v)
    {
        cout << elem << " ";
    }
    cout << endl;
}

// Faz a alocação de uma matriz
template <typename T>
vector<vector<T>> alocaMatriz(int m, int n)
{
    return vector<vector<T>>(m, vector<T>(n));
}

// Faz a alocação de uma matriz inicializada com um valor especifico
template <typename T>
vector<vector<T>> alocaMatriz(int m, int n, T x)
{
    return vector<vector<T>>(m, vector<T>(n, x));
}

// Faz a alocação de um vetor
template <typename T>
vector<T> alocaVetor(int n)
{
    return vector<T>(n);
}

// Faz a leitura do arquivo de entrada e o coloca em uma string com as linhas dividadas por \n
string lerTxt()
{
    ifstream file("input.txt");
    if (!file)
    {
        cout << "num abriu :(" << endl;
    }
    string str, input;
    while (getline(file, str))
    {

        input += str;
        input.push_back('\n');
    }
    return input;
}

// Calcula o número de variáveis sem contar as de folga
int calcularNumVariaveis(string input)
{
    int i = 0, num = 0;
    while (input[i] != '=')
    {
        i++;
    }
    while (input[i] != '\n')
    {
        if (input[i] == 'x')
        {
            num++;
        }
        i++;
    }
    return num;
}

/**
 * Calcula o tamanho da matriz A com base no número de variáveis e restrições.
 * @param input     [IN]  A string que possuí o input de onde os coeficientes serão lidos.
 * @param i    [OUT]  O número de variaveis contando as de folg.
 * @param j    [OUT]  O número de restrições.
 */
void calcularIJ(string input, int *i, int *j)
{
    int k = 0;
    (*i) = calcularNumVariaveis(input);
    (*j) = 0;

    while (input[k] != '\n')
    {
        k++;
    }

    while (input[k] != '\0')
    {
        if (input[k] == '>' || input[k] == '<')
        {
            (*i)++;
            (*j)++;
            if (input[k + 1] == '=')
            {
                k++;
            }
        }
        else if (input[k] == '=')
        {
            (*j)++;
        }
        k++;
    }
}

/**
 * Lê os coeficientes das restrições e as variáveis de folga.
 * @param A     [OUT]  A matriz contendo os coeficientes e as variáveis de folga.
 * @param input     [IN]  A string que possuí o input de onde os coeficientes serão lidos.
 * @param numRestricoes    [IN]  O número de restrições presentes no input.
 * @param numVariaveis    [IN]  O número total de variáveis presentes no input, contando as de folga.
 * @param b    [OUT]  O vetor de recursos.
 * @param c    [OUT]  O vetor de custos.
 */
void lerCoeficientes(vector<vector<double>> &A, string input, int numRestricoes, int numVariaveis, vector<double> &b, vector<double> &c)
{
    regex termoRegex("(-?\\d*[.]?\\d*)x(\\d+)");
    regex restricaoRegex("(<=|>=|<|>)");
    regex valRestricao(("\\d+[.]?\\d*"));
    smatch match;
    int linha = 0, k = 0, varFolga = 0;

    string custos;
    while (input[k] != '\n' && input[k] != '\0')
    {
        custos.push_back(input[k]);
        k++;
    }
    k++;

    for (int i = 0; i < numVariaveis; i++)
    {
        c[i] = 0;
    }

    auto it = custos.cbegin();
    while (regex_search(it, custos.cend(), match, termoRegex))
    {
        double coef = (match.str(1).empty() ? 1 : (match.str(1) == "-" ? -1 : stod(match.str(1))));
        int varIndex = stoi(match.str(2)) - 1;
        c[varIndex] = coef;
        it = match.suffix().first;
    }

    // Inicializa com zeros
    for (int i = 0; i < numRestricoes; i++)
    {
        for (int j = 0; j < numVariaveis; j++)
        {
            A[i][j] = 0;
        }
    }

    while (linha < numRestricoes)
    {
        // Separa a restrição
        string restricao;
        while (input[k] != '\n' && input[k] != '\0')
        {
            restricao.push_back(input[k]);
            k++;
        }
        k++; // Pular o '\n'

        auto it = restricao.cbegin();
        while (regex_search(it, restricao.cend(), match, termoRegex))
        {
            double coef = (match.str(1).empty() ? 1 : (match.str(1) == "-" ? -1 : stod(match.str(1))));
            int varIndex = stoi(match.str(2)) - 1;
            A[linha][varIndex] = coef;
            it = match.suffix().first;
        }

        // Separa os valores de b da restrição
        if (regex_search(it, restricao.cend(), match, valRestricao))
        {
            b[linha] = stod(match.str());
        }

        // Detecta o tipo da inequação e define variável de folga
        if (regex_search(restricao, match, restricaoRegex))
        {
            string sinal = match.str(1);
            if (sinal == "<=" || sinal == "<")
            {
                A[linha][calcularNumVariaveis(input) + varFolga] = 1.0;
            }
            else if (sinal == ">=" || sinal == ">")
            {
                A[linha][calcularNumVariaveis(input) + varFolga] = -1.0;
            } // Para '=' não adiciona nada (já está inicializado com zero)
            varFolga++;
        }

        linha++;
    }
}

/**
 * Separa uma submatriz de tamanho n-1 x n-1 da matriz M.
 * @param M     [IN]  Matriz original.
 * @param n     [IN]  Tamanho da matriz original.
 * @param i     [IN]  Linhas para serem removidas.
 * @param j     [IN]  Colunas para serem removidas.
 * @retval **R  Matriz reduzida resultante.
 */
template <typename T>
vector<vector<T>> matrizParcial(vector<vector<T>> M, int n, int m, vector<int> i, vector<int> j)
{
    vector<vector<T>> R = alocaMatriz<T>(n - 1, n - 1);

    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        if (!contemValor(i, k))
        {
            int coluna = 0;
            for (int q = 0; q < m; q++)
            {
                if (!contemValor(j, q))
                {
                    R[linha][coluna] = M[k][q];
                    coluna++;
                }
            }
            linha++;
        }
    }
    return R;
}

/**
 * Separa uma submatriz de tamanho n-1 x n-1 da matriz M, removendo apenas as colunas.
 * @param M     [IN]  Matriz original.
 * @param n     [IN]  Tamanho da matriz original.
 * @param j     [IN]  Colunas para serem removidas.
 * @retval **R  Matriz reduzida resultante.
 */
template <typename T>
vector<vector<T>> matrizParcialColuna(vector<vector<T>> M, int n, int m, vector<int> j)
{
    cout << "biava0.1" << endl;
    vector<vector<T>> R = alocaMatriz<T>(n - 1, n - 1);
    cout << "biava0.2" << endl;
    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        int coluna = 0;
        for (int q = 0; q < m; q++)
        {
            if (!contemValor(j, q))
            {
                R[linha][coluna] = M[k][q];
                coluna++;
                cout << "biava0.3" << endl;
            }
        }
        linha++;
    }
    cout << "biava0.4" << endl;
    return R;
}

/**
 * Separa uma submatriz de tamanho n-1 x n-1 da matriz M, removendo apenas as linhas.
 * @param M     [IN]  Matriz original.
 * @param n     [IN]  Tamanho da matriz original.
 * @param i     [IN]  Linhas para serem removidas.
 * @retval **R  Matriz reduzida resultante.
 */
template <typename T>
vector<vector<T>> matrizParcialLinha(vector<vector<T>> M, int n, int m, vector<int> i)
{
    vector<vector<T>> R = alocaMatriz<T>(n - 1, n - 1);

    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        if (!contemValor(i, k))
        {
            int coluna = 0;
            for (int q = 0; q < n; q++)
            {
                R[linha][coluna] = M[k][q];
                coluna++;
            }
            linha++;
        }
    }
    return R;
}

/**
 * Realiza o cálculo do determinante de uma matriz quadrada.
 * @param M     [IN]  Matriz que será usada para o cálculo.
 * @param n     [IN]  Tamanho da matriz.
 * @retval double  Resultado do determinante.
 */
double determinante(vector<vector<double>> M, int n)
{
    if (n == 1)
    {
        return M[0][0];
    }
    else
    {
        double result = 0;
        vector<vector<double>> matrizP;
        int i = 0;

        for (int j = 0; j < n; j++)
        {
            vector<int> linhas, colunas;
            linhas.push_back(i);
            colunas.push_back(j);
            matrizP = matrizParcial(M, n, linhas, colunas);
            result += pow(-1, i + j) * M[i][j] * determinante(matrizP, n - 1);
        }

        return result;
    }
}

/**
 * Realiza a multiplicação de duas matrizes.
 * @param A     [IN]  Matriz A.
 * @param lA    [IN]  Número de linhas da matriz A.
 * @param cA    [IN]  Número de colunas da matriz A.
 * @param B     [IN]  Matriz B.
 * @param cB    [IN]  Número de colunas da matriz B.
 * @retval **R  Matriz resultante da multiplicação.
 */
vector<vector<double>> multMatriz(vector<vector<double>> A, int lA, int cA, vector<vector<double>> B, int cB)
{
    vector<vector<double>> matrizFinal;
    matrizFinal = alocaMatriz<double>(lA, cB);
    for (int i = 0; i < lA; i++)
    {
        for (int j = 0; j < cB; j++)
        {
            double somaMult = 0;
            for (int k = 0; k < cA; k++)
            {
                somaMult = somaMult + (A[i][k] * B[k][j]);
            }
            matrizFinal[i][j] = somaMult;
        }
    }
    return matrizFinal;
}

/**
 * Realiza a multiplicação de uma matriz por um vetor.
 * @param A     [IN]  Matriz A.
 * @param lA    [IN]  Número de linhas da matriz A.
 * @param cA    [IN]  Número de colunas da matriz A.
 * @param B     [IN]  Vetor B.
 * @retval **R  Vetor resultante da multiplicação.
 */
vector<double> multMatrizVetor(vector<vector<double>> A, int lA, int cA, vector<double> B)
{
    vector<double> vetorFinal = alocaVetor<double>(lA);
    for (int i = 0; i < lA; i++)
    {
        double somaMult = 0;
        for (int k = 0; k < cA; k++)
        {
            somaMult = somaMult + (A[i][k] * B[k]);
        }
        vetorFinal[i] = somaMult;
    }
    return vetorFinal;
}

/**
 * Realiza a multiplicação de um vetor por uma matriz.
 * @param A     [IN]  Vetor A.
 * @param cA    [IN]  Número de colunas do vetor A.
 * @param B     [IN]  Matriz B.
 * @param cB    [IN]  Número de colunas da matriz B.
 * @retval **R  Vetor resultante da multiplicação.
 */
vector<double> multVetorMatriz(vector<double> A, int cA, vector<vector<double>> B, int cB)
{
    vector<double> vetorFinal = alocaVetor<double>(cA);

    for (int j = 0; j < cB; j++)
    {
        double somaMult = 0;
        for (int k = 0; k < cA; k++)
        {
            somaMult = somaMult + (A[k] * B[k][j]);
        }
        vetorFinal[j] = somaMult;
    }
    return vetorFinal;
}

/**
 * Realiza a multiplicação de dois vetores.
 * @param A     [IN]  Vetor A.
 * @param B     [IN]  Vetor B.
 * @retval double  Resultado da multiplicação.
 */
double multVetor(vector<double> A, vector<double> B)
{
    double somaMult = 0;
    for (int i = 0; i < A.size(); i++)
    {
        somaMult += A[i] * B[i];
    }
    return somaMult;
}

/**
 * Calcula a inversa de uma matriz quadrada.
 * @param M     [IN]  Matriz a ser invertida.
 * @param n     [IN]  Tamanho da matriz.
 * @retval **R  Matriz inversa resultante.
 */
vector<vector<double>> inversa(vector<vector<double>> M, int n)
{
    cout << "biava0.5" << endl;
    vector<vector<double>> I;
    cout << "biava" << endl;
    I = alocaMatriz<double>(n, n, 0.0);
    cout << "biava2" << endl;
    for (int i = 0; i < n; i++)
    {
        I[i][i] = 1.0;
    }

    for (int i = 0; i < n; i++)
    {
        double divisor = M[i][i];
        if (divisor == 0.0)
        {
            throw runtime_error("Matriz singular: divisão por zero.");
        }
        for (int j = 0; j < n; j++)
        {
            M[i][j] /= divisor;
            I[i][j] /= divisor;
        }

        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                double multiplier = M[j][i];
                for (int k = 0; k < n; k++)
                {
                    M[j][k] -= multiplier * M[i][k];
                    I[j][k] -= multiplier * I[i][k];
                }
            }
        }
    }
    return I;
}

/**
 * Escolhe colunas aleatórias para formar a matriz básica.
 * @param M     [IN]  Matriz original.
 * @param m     [IN]  Número de linhas da matriz.
 * @param n     [IN]  Número de colunas da matriz.
 * @param B     [OUT]  Vetor que armazenará as colunas da matriz basica.
 * @param N     [OUT]  Vetor que armazenará as colunas da matriz não basica.
 */
void escolherColunasAleatorias(vector<vector<double>> M, int m, int n, vector<int> &B, vector<int> &N)
{
    vector<vector<double>> MatrizBasicaTemp = alocaMatriz<double>(m, m, 0.0);
    vector<vector<int>> conjuntosTestados = alocaMatriz<int>(n, m);
    int ValBTemp;

    while (determinante(MatrizBasicaTemp, m) == 0)
    {
        vector<int> BTemp(m, -1);
        for (int i = 0; i < m; i++)
        {
            bool existe = true;
            cout << "Escolhendo coluna " << i << endl;
            while (existe)
            {
                ValBTemp = (rand() % n);
                if (contemValor(BTemp, ValBTemp))
                {
                    existe = true;
                }
                else
                {
                    existe = false;
                }
            }
            BTemp[i] = ValBTemp;
            B[i] = ValBTemp;
            for (int j = 0; j < m; j++)
            {
                cout << M[j][B[i]] << endl;
                MatrizBasicaTemp[j][i] = M[j][B[i]];
            }
        }
        if (!contemConjunto(conjuntosTestados, B))
        {
            conjuntosTestados.push_back(B);
            cout << "Matriz basica: " << endl;
            imprimeMatriz(MatrizBasicaTemp);
            cout << "det: " << determinante(MatrizBasicaTemp, m) << endl;
        }
    }
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        cout << i << endl;
        if (!contemValor(B, i))
        {
            N[j] = i;
            j++;
        }
    }
}


vector<double> calcSolucaoBasica(vector<vector<double>> A, vector<int> B, vector<int> N, vector<double> b)
{
    vector<double> solucaoBasica;
    vector<vector<double>> inversaB = inversa(matrizParcialColuna(A, A.size(), N), B.size());
    solucaoBasica = multMatrizVetor(inversaB, B.size(), B.size(), b);
    return solucaoBasica;
}

double calcCustosRelativos(vector<vector<double>> A, vector<int> B, vector<int> N, vector<double> c)
{
    vector<double> vetorMultiplicador, custosRelativos, custosBasica;
    vector<vector<double>> inversaB = inversa(matrizParcialColuna(A, A.size(), N), B.size());
    custosBasica = alocaVetor<double>(B.size());
    for (int i = 0; i < B.size(); i++)
    {
        custosBasica[i] = c[B[i]];
    }
    vetorMultiplicador = multVetorMatriz(custosBasica, custosBasica.size(), inversaB, B.size());
    custosRelativos = alocaVetor<double>(N.size());
    for(int i = 0; i < N.size(); i++)
    {
        custosRelativos[i] = c[N[i]] - multVetor(vetorMultiplicador, A[N[i]]);
    }
    return menorValor(custosRelativos);
}

/**
 * Realiza a fase II do método Simplex.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazenará as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazenará as colunas da matriz não basica.
 * @param b     [IN]  Vetor de recursos.
 * @param c     [IN]  Vetor de custos.
 */
void faseII(vector<vector<double>> A, vector<int> &B, vector<int> &N, vector<double> b, vector<double> c)
{
    vector<double> solucaoBasica, y;
    cout << "pre-inversa" << endl;
    imprimeMatriz(matrizParcialColuna(A, A.size(), N));
    vector<vector<double>> inversaB = inversa(matrizParcialColuna(A, A.size(), N), B.size());
    double custoRelativo = -1;
    int k = 0;
    cout << "Antes" << endl;
    while (custoRelativo < 0 && k < N.size())
    {
        cout << "k: " << k << endl;
        cout << "B: " << endl;
        impremeVetor(B);
        cout << "N: " << endl;
        impremeVetor(N);
        solucaoBasica = calcSolucaoBasica(A, B, N, b);
        custoRelativo = calcCustosRelativos(A, B, N, c);
        y = multMatrizVetor(inversaB, B.size(), B.size(), A[N[k]]);
        if(vetorMenorZero(y)){
            throw runtime_error("Solução não é viável.");
            return;
        }
        int variavelSaida;
        vector<double> aux = alocaVetor<double>(y.size());
        for (int i = 0; i < y.size(); i++)
        {
            if(y[i] != 0){
                aux[i] = solucaoBasica[i] / y[i];
            }
            else {
                aux[i] = INT16_MAX;
            }
        }
        double e = menorValor(aux);
        for (int i = 0; i < y.size(); i++)
        {
            if (aux[i] == e)
            {
                variavelSaida = i;
                break;
            }
        }
        int varAux = B[variavelSaida];
        B[variavelSaida] = N[k];
        N[k] = varAux;
        k++;
    }
}

int main()
{
    string input = lerTxt();
    int numVariaveis, numRestricoes;
    calcularIJ(input, &numVariaveis, &numRestricoes);

    cout << input << endl;

    vector<vector<double>> A = alocaMatriz<double>(numRestricoes, numVariaveis);
    vector<double> b(numRestricoes), c(numVariaveis);

    lerCoeficientes(A, input, numRestricoes, numVariaveis, b, c);

    cout << "Vetor c: " << endl;
    impremeVetor(c);

    cout << endl
         << "Matriz A: " << endl;
    imprimeMatriz(A);

    cout << endl
         << "Vetor b: " << endl;
    impremeVetor(b);
    cout << endl;

    vector<int> B(numRestricoes), N(numVariaveis - numRestricoes);
    escolherColunasAleatorias(A, numRestricoes, numVariaveis, B, N);
    cout << "escolheu as aleatorias" << endl;
    faseII(A, B, N, b, c);

    cout << endl
         << "B: " << endl;
    impremeVetor(B);
    cout << endl;
    cout << "N: " << endl;
    impremeVetor(N);
    cout << endl;
    /*
    vector<vector<double>> M = alocaMatriz<double>(3, 3);
    M[0][0] = 2;
    M[0][1] = 1;
    M[0][2] = 4;
    M[1][0] = 0;
    M[1][1] = 2;
    M[1][2] = 1;
    M[2][0] = 3;
    M[2][1] = 0;
    M[2][2] = 5;
    cout << determinante(M, 3) << endl;

    vector<vector<double>> Ml;

    Ml = inversa(M, 3);
    imprimeMatriz(M);

    cout << endl << endl << "inversa: " << endl;

    imprimeMatriz(Ml);

    cout << "Jean" << endl << endl << endl;
    vector<vector<double>> I;
    cout << "Jean1.5" << endl << endl << endl;
    I = alocaMatriz<double>(3, 3, 0.0);
    cout << "Jean2" << endl << endl << endl;
    I[0][0] = 1.0;
    I[1][1] = 1.0;
    I[2][2] = 1.0;

    cout << "biava" << endl;
    M = multMatriz(M, 3, 3, I, 3);

    cout << "tashuka" << endl;
    imprimeMatriz(M);



    vector<vector<double>> N = alocaMatriz<double>(2, 2);
    N[0][0] = 1;
    N[0][1] = 5;
    N[1][0] = 2;
    N[1][1] = 3;
    cout << determinante(N, 2) << endl;
    vector<vector<double>> F = alocaMatriz<double>(2, 2);
    cout << "foi" << endl;
    */
    return 0;
}
