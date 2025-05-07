#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <math.h>
#include <map>
#include <vector>

using namespace std;

// Faz a alocação de uma matriz
template <typename T>
vector<vector<T>>  alocaMatriz(int m, int n)
{
    try {
        cout << "Alocando matriz de " << m << " x " << n << endl;
        return vector<vector<T>>(m, vector<T>(n));  
    } catch (const std::bad_alloc& e) {
        std::cerr << "Erro de alocação: " << e.what() << std::endl;
        return vector<vector<T>>(m, vector<T>(n));
    }
    
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

    for(int i = 0; i < numVariaveis; i++){
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
 * @param i     [IN]  Linha para ser removida.
 * @param j     [IN]  Coluna para ser removida.
 * @retval **R  Matriz reduzida resultante.
 */
template <typename T>
vector<vector<T>> matrizParcial(vector<vector<T>> M, int n, int i, int j)
{
    vector<vector<T>> R = alocaMatriz<T>(n - 1, n - 1);

    int linha = 0;
    for (int k = 0; k < n; k++)
    {
        if (k == i)
        {
            continue;
        }
        int coluna = 0;
        for (int q = 0; q < n; q++)
        {
            if (q == j)
            {
                continue;
            }
            R[linha][coluna] = M[k][q];
            coluna++;
        }
        linha++;
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
            matrizP = matrizParcial(M, n, i, j);
            result += pow(-1, i + j) * M[i][j] * determinante(matrizP, n - 1);
        }
        
        return result;
    }
}

vector<vector<double>> multMatriz(vector<vector<double>> A, int lA, int cA, vector<vector<double>> B, int cB){
    vector<vector<double>> matrizFinal;
    matrizFinal = alocaMatriz<double>(lA, cB);
    for(int i=0; i<lA; i++){
        for(int j=0; j<cB; j++){
            double somaMult=0;
            for(int k=0; k<cA; k++){
                somaMult = somaMult + (A[i][k] * B[k][j]);
            }
            matrizFinal[i][j] = somaMult;
        }
    }
    return matrizFinal;
}

vector<vector<double>> inversa(vector<vector<double>> M, int n){
    vector<vector<double>> I;
    I = alocaMatriz<double>(n, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i == j){
                I[i][j] = 1;
            }
            else {
                I[i][j] = 0;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        double divisor = M[i][i];
        for (int j = 0; j <= n; j++){
            M[i][j] /= divisor;
            I[i][j] /= divisor;
        }

        for (int j = 0; j < n; j++){
            if (j != i) {
                double multiplier = M[j][i];
                for (int k = 0; k <= n; k++){
                    M[j][k] -= multiplier * M[i][k];
                    I[j][k] -= multiplier * I[i][k];
                }      
            }
        }
    }
    return I;

}

int fatorial(int n){
    if (n == 1){
        return n;
    }
    else {
        return n * fatorial(n - 1);
    }
}

/*
void escolherColunas(double **M, int m, int n, int *B, int *N){
    int conjuntosTestados[fatorial(n)][m], numConjuntos, colunasTemp[m];
    double **MatrizTemp;
    for (int i = 0; i < m; i++){
        colunasTemp[i] = i;
    }
    for (int coluna : colunasTemp){

    }

        B[0] = ((m - 1) % rand()) + 1;
        int aux = B[0];
        bool existe = true;
        for (int i = 1; i < m; i++){
            while (existe){
                aux = ((m - 1) % rand()) + 1;
                int j = 0;
                while (j < i && !existe){
                    if (B[j] == aux){
                        existe = true; 
                    }
                    else {
                        existe = false;
                    }
                    j++;
                }
            }
            B[i] = aux;
        }
    
}
*/

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

    cout << endl << endl;

    imprimeMatriz(Ml);

    cout << "Jean" << endl << endl << endl;
    
    vector<vector<double>> I;
    cout << "Jean1.5" << endl << endl << endl;
    I = alocaMatriz<double>(3, 3);
    cout << "Jean2" << endl << endl << endl;
    I[0][0] = 1;
    I[0][1] = 0;
    I[0][2] = 0;
    I[1][0] = 0;
    I[1][1] = 1;
    I[1][2] = 0;
    I[2][0] = 0;
    I[2][1] = 0;
    I[2][2] = 1;
    
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

    return 0;
}

