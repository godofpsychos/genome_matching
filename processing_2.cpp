#include <bits/stdc++.h>
using namespace std;
enum Column
{
    Start_Index, // Start_Index can be represented by 1
    Read,
    Fred
};
enum ReadType
{
    Match,
    Unmatch
};
std::string getSecondLine(const std::string &filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "Could not open the file!" << std::endl;
        return "";
    }

    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line))
    {
        lineNumber++;
        if (lineNumber == 2)
        {
            file.close();
            return line; // Return the second line
        }
    }

    file.close();
    return ""; // Return an empty string if the second line doesn't exist
}
double logSumExp(const std::vector<double> &x)
{
    double max_val = *std::max_element(x.begin(), x.end());
    double sum = 0.0;

    for (double val : x)
    {
        sum += exp(val - max_val);
    }

    return max_val + log(sum);
}
double calculate_pval(vector<double> err_arr, int N, int K)
{
    vector<double> pl(K), pl_prev(K);
    pl_prev[0] = 0;
    double ln_pval, ln_pval_prev;
    for (int n = 1; n <= N; n++)
    {
        double pn = err_arr[n - 1];
        double ln_pn = log(pn);
        double ln_1_pn = log(1 - pn);
        if (n < K)
            pl_prev[n] = -1e100;
        int bound = n < K - 1 ? n : K - 1;
        for (int k = 1; k <= bound; k++)
        {
            pl[k] = logSumExp({pl_prev[k - 1] + ln_pn, pl_prev[k] + ln_1_pn});
        }
        pl[0] = pl_prev[0] + ln_1_pn;
        if (n == K)
            ln_pval = pl_prev[K - 1] + ln_pn;
        else if (n > K)
            ln_pval = logSumExp({pl_prev[K - 1] + ln_pn, ln_pval_prev});
        pl_prev = pl;
        ln_pval_prev = ln_pval;
    }
    return ln_pval;
}

int main()
{
    const string ref_gnome_file_path = "/home/tarun/Desktop/genome_matching/EPI_ISL_402124-ORF1ab.fasta",
                 processed_data_csv_file_path = "/home/tarun/Desktop/genome_matching/processed_bam_data.csv";
    string ref_gnome = getSecondLine(ref_gnome_file_path);
    int len_ref_genome = ref_gnome.size();
    vector<vector<int>> processedReadData(len_ref_genome, vector<int>(2, 0));
    vector<vector<double>> colmwise_fred(len_ref_genome);
    std::ifstream file(processed_data_csv_file_path);
    vector<std::string> row;
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Could not open the file!" << std::endl;
        return 0;
    }
    if (!std::getline(file, line))
    {
        throw("Given Processed BAM File is Empty...!");
    }
    for (int col = 0; col < len_ref_genome; col++)
    {
        while (std::getline(file, line))
        {
            row.clear();
            std::stringstream ss(line);
            std::string cell;

            // Split the line into columns
            while (std::getline(ss, cell, '|'))
            {
                row.push_back(cell);
            }

            // Process the row if it has sufficient columns
            if (row.size() > 1)
            {
                // cout<<row[Column::Read]<<"\n";
                int start_idx = std::stoi(row[Column::Start_Index], nullptr, 10);
                const std::string &read = row[Column::Read];
                vector<double> fred_values;
                std::string fred_column = row[Column::Fred]; // The "Fred" column contains comma-separated values

                std::stringstream fred_stream(fred_column);
                std::string fred_value;

                // Split the "Fred" column values by commas and convert to double
                while (std::getline(fred_stream, fred_value, ','))
                {
                    cout<<fred_value<<"\n";
                    // double value = std::stod(fred_value); // Convert the value to double
                    // fred_values.push_back(value);         // Store it in the vector
                }

                // Only process if the read is non-empty and the start index is valid
                if (start_idx >= 0 && !read.empty())
                {
                    int len_read = read.size();

                    // Pre-allocate the size of processedReadData (if not already done) for efficiency
                    for (int i = 0; i < len_read; ++i)
                    {
                        int position = start_idx + i;
                        // Ensure the vector size is large enough
                        if (position >= processedReadData.size())
                        {
                            processedReadData.resize(position + 1, std::vector<int>(2, 0)); // Resize with 2 counters: match & unmatch
                        }

                        if (ref_gnome[position] == read[i])
                        {
                            processedReadData[position][ReadType::Match]++; // match++
                        }
                        else
                        {
                            processedReadData[position][ReadType::Unmatch]++; // unmatch++
                        }
                        colmwise_fred[position].push_back(fred_values[i]);
                    }
                }
            }
        }
    }

    file.close();

    cout << "Data Processed Sucessfully\n";

    return 0;
}