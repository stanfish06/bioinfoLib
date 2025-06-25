#include <vector>
using namespace std;
// #include <iostream>

// this function checks if Ax = b is solvable with mod2 arithmetics
bool mod2Solve(vector<vector<bool>> Ab)
{
    int nr = Ab.size();
    int nc = Ab[0].size();
    int row = 0;
    for (int i = 0; i < (nc - 1); ++i)
    {
        int pivot = -1;
        for (int j = row; j < nr; ++j)
        {
            if (Ab[j][i])
            {
                pivot = j;
                break;
            }
        }

        if (pivot != -1)
        {
            if (pivot != row)
            {
                swap(Ab[row], Ab[pivot]);
            }
            for (int j = 0; j < nr; ++j)
            {
                if (j != row && Ab[j][i])
                {
                    for (int k = 0; k < nc; k++)
                    {
                        Ab[j][k] = Ab[j][k] ^ Ab[row][k];
                    }
                }
            }
            ++row;
        }
    }

    // rows with pivots are defintiely solvable, so we check others.
    for (int i = (row - 1); i < nr; ++i)
    {
        if (Ab[i][nc - 1])
        {
            bool allzero = true;
            for (int j = 0; j < (nc - 1); ++j)
            {
                if (Ab[i][j])
                {
                    allzero = false;
                    break;
                }
            }
            if (allzero)
            {
                return false;
            }
        }
    }

    return true;
}

int main()
{
    // cout << "demo:" << endl;

    // // --- Test case 1 ---
    // // 1 1 0 | 1
    // // 1 0 1 | 1
    // // 0 1 1 | 0
    // // Expected: Solvable
    // int nr1 = 3, nc1 = 3;
    // vector<vector<bool>> Ab1(nr1, vector<bool>(nc1 + 1, false));
    // Ab1[0][0] = true;
    // Ab1[0][1] = true;
    // Ab1[0][nc1] = true; // Row 0: 1 1 0 | 1
    // Ab1[1][0] = true;
    // Ab1[1][2] = true;
    // Ab1[1][nc1] = true; // Row 1: 1 0 1 | 1
    // Ab1[2][1] = true;
    // Ab1[2][2] = true; // Row 2: 0 1 1 | 0

    // cout << "System 1:" << endl;
    // for (int i = 0; i < nr1; ++i)
    // {
    //     for (int j = 0; j < nc1; ++j)
    //     {
    //         cout << Ab1[i][j] << " ";
    //     }
    //     cout << "| " << Ab1[i][nc1] << endl;
    // }
    // // Call your existing mod2Solve function.
    // // Pass a copy if mod2Solve modifies the input vector.
    // vector<vector<bool>> Ab1_copy = Ab1;
    // cout << "Test case 1: " << (mod2Solve(Ab1_copy) ? "Solvable" : "Not Solvable") << endl
    //      << endl;

    // // --- Test case 2 ---
    // // 1 0 0 | 1
    // // 1 1 1 | 0
    // // 0 1 1 | 0
    // // Expected: Not Solvable
    // int nr2 = 3, nc2 = 3;
    // vector<vector<bool>> Ab2(nr2, vector<bool>(nc2 + 1, false));
    // Ab2[0][0] = true;
    // Ab2[0][nc2] = true; // Row 0: 1 0 0 | 1
    // Ab2[1][0] = true;
    // Ab2[1][1] = true;
    // Ab2[1][2] = true; // Row 1: 1 1 1 | 0
    // Ab2[2][1] = true;
    // Ab2[2][2] = true;
    // // Ab2[2][nc2] = true; // Row 2: 0 1 1 | 1

    // cout << "System 2:" << endl;
    // for (int i = 0; i < nr2; ++i)
    // {
    //     for (int j = 0; j < nc2; ++j)
    //     {
    //         cout << Ab2[i][j] << " ";
    //     }
    //     cout << "| " << Ab2[i][nc2] << endl;
    // }
    // // Call your existing mod2Solve function.
    // // Pass a copy if mod2Solve modifies the input vector.
    // vector<vector<bool>> Ab2_copy = Ab2;
    // cout << "Test case 2: " << (mod2Solve(Ab2_copy) ? "Solvable" : "Not Solvable") << endl;

    return 0;
}