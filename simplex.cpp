#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <math.h>
#include <map>
#include <vector>
#include <chrono>

using namespace std;

/**
 * Extrai uma coluna específica de uma matriz.
 * @param A     [IN]  A matriz de onde a coluna será extraída.
 * @param i     [IN]  O índice da coluna a ser extraída.
 * @retval vector<double>  Retorna um vetor contendo os elementos da coluna.
 */
vector<double> pegaColuna(const vector<vector<double>> &A, int i)
{
    vector<double> coluna(A.size());
    for (int j = 0; j < A.size(); ++j)
    {
        coluna[j] = A[j][i];
    }
    return coluna;
}

/**
 * Extrai colunas específicas de uma matriz.
 * @param A     [IN]  A matriz de onde as colunas serão extraídas.
 * @param B     [IN]  Um vetor contendo os índices das colunas a serem extraídas.
 * @retval vector<vector<double>>  Retorna uma matriz contendo as colunas extraídas.
 */
vector<vector<double>> extraiColunas(const vector<vector<double>> A, const vector<int> B)
{
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

// Faz a impressão de uma matriz
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

// Faz a impressão de um vetor
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

// Faz a alocação de um vetor
template <typename T>
vector<T> alocaVetor(int n, T x)
{
    return vector<T>(n, x);
}

// Faz a leitura do arquivo de entrada e o coloca em uma string com as linhas dividadas por \n
string lerTxt()
{
    ifstream file("input.txt");
    if (!file)
    {
        throw runtime_error("Não foi possivel abrir o arquivo.");
    }
    string str, input;
    while (getline(file, str))
    {
        input += str;
        input.push_back('\n');
    }
    for (int i = 0; i < input.size(); i++)
    {
        if (input[i] == ' ')
        {
            input.erase(i, 1);
            i--;
        }
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
    regex termoRegex("(-?\\d*[.]?\\d*)x(\\d+)", regex::optimize);
    regex restricaoRegex("(<=|>=|<|>)", regex::optimize);
    regex valRestricao(("\\d+[.]?\\d*"), regex::optimize);
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
vector<vector<T>> matrizParcial(vector<vector<T>> M, int n, vector<int> i, vector<int> j)
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
vector<vector<T>> matrizParcialColuna(vector<vector<T>> M, int n, vector<int> j)
{
    vector<vector<T>> R = alocaMatriz<T>(n - 1, n - 1);

    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        int coluna = 0;
        for (int q = 0; q < n; q++)
        {
            if (!contemValor(j, q))
            {
                R[linha][coluna] = M[k][q];
                coluna++;
            }
        }
        linha++;
    }
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
vector<vector<T>> matrizParcialLinha(vector<vector<T>> M, int n, vector<int> i)
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
 * @retval **R  Matriz inversa resultante.
 */
vector<vector<double>> inversa(vector<vector<double>> M)
{
    int n = M.size();
    vector<vector<double>> I = alocaMatriz<double>(n, n, 0.0);

    // Inicializa a matriz identidade
    for (int i = 0; i < n; i++)
        I[i][i] = 1.0;

    for (int i = 0; i < n; i++)
    {
        // Procura o maior elemento na coluna para evitar divisão por zero
        double maxEl = fabs(M[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++)
        {
            if (fabs(M[k][i]) > maxEl)
            {
                maxEl = fabs(M[k][i]);
                maxRow = k;
            }
        }
        if (fabs(maxEl) == 0.0)
            throw runtime_error("Matriz singular, não possui inversa.");

        // Troca as linhas na matriz original e na identidade
        swap(M[i], M[maxRow]);
        swap(I[i], I[maxRow]);

        // Divide a linha pelo pivô
        double piv = M[i][i];
        for (int j = 0; j < n; j++)
        {
            M[i][j] /= piv;
            I[i][j] /= piv;
        }

        // Elimina os outros elementos da coluna
        for (int k = 0; k < n; k++)
        {
            if (k != i)
            {
                double f = M[k][i];
                for (int j = 0; j < n; j++)
                {
                    M[k][j] -= f * M[i][j];
                    I[k][j] -= f * I[i][j];
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
                MatrizBasicaTemp[j][i] = M[j][B[i]];
            }
        }
        if (!contemConjunto(conjuntosTestados, B))
        {
            conjuntosTestados.push_back(B);
        }
    }
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        if (!contemValor(B, i) && !contemValor(N, i))
        {
            N[j] = i;
            j++;
        }
    }
}

/**
 * Calcula a solução básica inicial.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazena as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazena as colunas da matriz não basica.
 * @param b     [IN]  Vetor de recursos.
 * @retval vector<double>  Retorna a solução básica inicial.
 */
vector<double> calcSolucaoBasica(vector<vector<double>> A, vector<int> B, vector<int> N, vector<double> b)
{
    vector<double> solucaoBasica = alocaVetor<double>(A[0].size(), 0.0);
    vector<vector<double>> inversaB = inversa(extraiColunas(A, B));
    vector<double> aux = multMatrizVetor(inversaB, B.size(), B.size(), b);
    for (int i = 0; i < B.size(); i++)
    {
        solucaoBasica[B[i]] = aux[i];
    }
    int j = 0;

    return solucaoBasica;
}

/**
 * Calcula os custos relativos e identifica a variável de entrada.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazena as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazena as colunas da matriz não basica.
 * @param c     [IN]  Vetor de custos.
 * @param varEntrada [OUT]  Variável de entrada identificada.
 * @retval double  Retorna o menor custo relativo encontrado.
 */
double calcCustosRelativos(vector<vector<double>> A, vector<int> B, vector<int> N, vector<double> c, int &varEntrada)
{
    vector<double> vetorMultiplicador, custosRelativos, custosBasica;
    vector<vector<double>> inversaB = inversa(extraiColunas(A, B));
    custosBasica = alocaVetor<double>(B.size());
    for (int i = 0; i < B.size(); i++)
    {
        custosBasica[i] = c[B[i]];
    }
    vetorMultiplicador = multVetorMatriz(custosBasica, custosBasica.size(), inversaB, B.size());
    custosRelativos = alocaVetor<double>(N.size());
    for (int i = 0; i < N.size(); i++)
    {
        custosRelativos[i] = c[N[i]] - multVetor(vetorMultiplicador, pegaColuna(A, N[i]));
    }
    double menor = menorValor(custosRelativos);
    for (int i = 0; i < custosRelativos.size(); i++)
    {
        if (custosRelativos[i] == menor)
        {
            varEntrada = i;
            break;
        }
    }
    return menor;
}

/**
 * Verifica se a fase I é necessária.
 * @param input     [IN]  A string que possuí o input de onde os coeficientes serão lidos.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazena as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazena as colunas da matriz não basica.
 * @param b     [IN]  Vetor de recursos.
 * @param c     [IN]  Vetor de custos.
 * @retval bool  Retorna true se a fase I for necessária, false caso contrário.
 */
bool preFaseI(string input, vector<vector<double>> A, vector<int> &B, vector<int> &N, vector<double> &b, vector<double> &c)
{
    regex max("max", regex_constants::icase, regex::optimize);
    smatch match;
    regex_search(input, match, max);
    if (!match.empty())
    {
        for (int i = 0; i < c.size(); i++)
        {
            c[i] = -c[i];
        }
    }
    for (int i = 0; i < b.size(); i++)
    {
        if (b[i] < 0)
        {
            b[i] = -b[i];
            for (int j = 0; j < A[i].size(); j++)
            {
                A[i][j] = -A[i][j];
            }
        }
    }
    int i = 0;
    while (input[i] != '\n' && input[i] != '\0')
    {
        i++;
    }
    regex restricaoRegex("(>=|>|=)", regex::optimize);
    auto it = input.cbegin();
    it += i + 1;
    while (regex_search(it, input.cend(), match, restricaoRegex))
    {
        string op = match.str(1);
        size_t pos = match.position(1) + (it - input.cbegin());
        if (op == "=")
        {
            // Verifica se não é parte de <= ou >=
            if (!(pos > 0 && (input[pos - 1] == '<' || input[pos - 1] == '>')))
            {
                return true; // '=' sozinho
            }
        }
        else if (op == ">" || op == ">=")
        {
            return true;
        }
        it += match.position(0) + match.length(0);
    }
    return false;
}

/**
 * Realiza a fase I do método Simplex.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazena as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazena as colunas da matriz não basica.
 * @param b     [IN]  Vetor de recursos.
 * @param c     [IN]  Vetor de custos.
 */

void faseI(vector<vector<double>> A, vector<int> &B, vector<int> &N, vector<double> b, vector<double> c)
{
    /*Fase 1 ta tudo errada
    cout << "Fase I (Simplex):" << endl;

    int m = b.size();
    int n = A[0].size();

    // Identifica restrições que precisam de variáveis artificiais
    vector<int> artificiais;
    vector<vector<double>> A_aux = A;
    vector<double> c_aux = alocaVetor<double>(n + m, 0.0);
    for (int i = n; i < n + m; ++i) {
        c_aux[i] = 1.0;
    }

    // Adiciona variáveis artificiais para cada restrição que não tem folga positiva
    for (int i = 0; i < m; ++i) {
        // Adiciona coluna artificial
        for(int j = 0; j < m; ++j) {
            if (A_aux[i][j] < 0) {
                A_aux[i].push_back(1.0); // Variável artificial
            } else {
                A_aux[i].push_back(0.0); // Coluna de zeros
            }
        }
    }

    N.resize(n);
    for(int i = 0; i < n; i++)
    {
        N[i] = i;
    }

    for(int i = 0; i < m; i++)
    {
        B[i] = i + n;
    }

    cout << "A_aux: " << endl;
    imprimeMatriz(A_aux);
    cout << "Básicas: "; impremeVetor(B);
    cout << "Não básicas: "; impremeVetor(N);
    cout << "c_aux: "; impremeVetor(c_aux);

    // Inicializa solução básica
    vector<double> solucaoBasica = calcSolucaoBasica(A_aux, B, N, b);

    cout << "A_aux: " << endl;
    imprimeMatriz(A_aux);
    cout << "Básicas: "; impremeVetor(B);
    cout << "Não básicas: "; impremeVetor(N);
    cout << "c_aux: "; impremeVetor(c_aux);

    // Simplex Fase I: Minimiza soma das artificiais
    int iter = 0;
    while (true) {
        int varEntrada;
        double custoRelativo = calcCustosRelativos(A_aux, B, N, c_aux, varEntrada);
        if (custoRelativo >= 0) break; // Ótimo encontrado

        // Calcula direção
        vector<vector<double>> invB = inversa(extraiColunas(A_aux, B));
        vector<double> y = multMatrizVetor(invB, B.size(), B.size(), pegaColuna(A_aux, N[varEntrada]));

        // Teste de viabilidade
        if (vetorMenorZero(y)) {
            throw runtime_error("Problema inviável na Fase I.");
        }

        // Razão mínima
        int varSaida = -1;
        double minRazao = std::numeric_limits<double>::max();
        for (int i = 0; i < y.size(); ++i) {
            if (y[i] > 1e-8) {
                double razao = solucaoBasica[B[i]] / y[i];
                if (razao < minRazao) {
                    minRazao = razao;
                    varSaida = i;
                }
            }
        }
        if (varSaida == -1) {
            throw runtime_error("Solução ilimitada na Fase I.");
        }

        // Troca básica
        int temp = B[varSaida];
        B[varSaida] = N[varEntrada];
        N[varEntrada] = temp;

        // Atualiza solução básica
        solucaoBasica = calcSolucaoBasica(A_aux, B, N, b);
        iter++;
    }

    // Verifica se todas as artificiais são zero
    for (int idx : artificiais) {
        if (fabs(solucaoBasica[idx]) > 1e-8) {
            throw runtime_error("Problema original inviável (artificiais não zeradas).");
        }
    }

    // Remove variáveis artificiais de B, N, A e c
    sort(artificiais.rbegin(), artificiais.rend());
    for (int idx : artificiais) {
        for (auto& row : A_aux) row.erase(row.begin() + idx);
        c_aux.erase(c_aux.begin() + idx);
        B.erase(remove(B.begin(), B.end(), idx), B.end());
        N.erase(remove(N.begin(), N.end(), idx), N.end());
    }

    cout << "Fase I concluída. Solução básica viável encontrada." << endl;
    cout << "Básicas: "; impremeVetor(B);
    cout << "Não básicas: "; impremeVetor(N);
    */
}

/**
 * Realiza a fase II do método Simplex.
 * @param A     [IN]  Matriz original.
 * @param B     [IN]  Vetor que armazena as colunas da matriz basica.
 * @param N     [IN]  Vetor que armazena as colunas da matriz não basica.
 * @param b     [IN]  Vetor de recursos.
 * @param c     [IN]  Vetor de custos.
 */
void faseII(vector<vector<double>> A, vector<int> &B, vector<int> &N, vector<double> b, vector<double> c)
{
    vector<double> solucaoBasica, y;
    double custoRelativo = -1;
    int k = 0;

    while (custoRelativo < 0)
    {
        vector<vector<double>> inversaB = inversa(extraiColunas(A, B));
        solucaoBasica = calcSolucaoBasica(A, B, N, b);
        int variavelEntrada;
        custoRelativo = calcCustosRelativos(A, B, N, c, variavelEntrada);
        if (custoRelativo < 0)
        {
            y = multMatrizVetor(inversaB, B.size(), B.size(), pegaColuna(A, N[variavelEntrada]));
            if (vetorMenorZero(y))
            {
                throw runtime_error("Solução não é viável.");
                return;
            }
            int variavelSaida;
            int menor = INT16_MAX;
            vector<double> aux = alocaVetor<double>(y.size());
            for (int i = 0; i < y.size(); i++)
            {
                if (y[i] > 0)
                {
                    aux[i] = solucaoBasica[B[i]] / y[i];
                    if (aux[i] < menor)
                    {
                        menor = aux[i];
                        variavelSaida = i;
                    }
                }
            }
            int varAux = B[variavelSaida];
            B[variavelSaida] = N[variavelEntrada];
            N[variavelEntrada] = varAux;
            k++;
            /*
            cout << "Custo relativo: " << custoRelativo << endl;
            cout << "Básica: " << endl;
            impremeVetor(B);
            cout << "Não básica: " << endl;
            impremeVetor(N);
            */
        }
    }

    cout << "Solucao otima encontrada." << endl;
    double solucaoOtima = 0;
    for (int i = 0; i < (B.size() + N.size()); i++)
    {
        solucaoOtima += c[i] * solucaoBasica[i];
        cout << "c[" << i << "] = " << c[i] << ", solucaoBasica[" << i << "] = " << solucaoBasica[i] << "  ";
    }
    cout << endl
         << "Solucao otima: " << solucaoOtima << endl;
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
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
    cout << endl
         << "-------------------------" << endl;

    vector<int> B(numRestricoes), N(numVariaveis - numRestricoes);

    if (preFaseI(input, A, B, N, b, c))
    {
        cout << "Fase I necessaria." << endl;
        faseI(A, B, N, b, c);
    }
    else
    {
        cout << "Fase I nao necessaria." << endl;
        escolherColunasAleatorias(A, numRestricoes, numVariaveis, B, N);
    }
    faseII(A, B, N, b, c);

    /*
    preFaseI(input, A, B, N, b, c);
    escolherColunasAleatorias(A, numRestricoes, numVariaveis, B, N);
    faseII(A, B, N, b, c);
    */
    cout << endl
         << "B: " << endl;
    impremeVetor(B);
    cout << endl;
    cout << "N: " << endl;
    impremeVetor(N);
    cout << endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Tempo demorado: " << duration.count() << " microsegundos" << std::endl;
    return 0;
}
