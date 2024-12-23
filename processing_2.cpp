#include <bits/stdc++.h>
using namespace std;
enum Column
{
    Start_Index, // Start_Index can be represented by 1
    Read
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
std::vector<std::string> getCSVHeader(const std::string &filePath)
{
    std::ifstream file(filePath);
    std::vector<std::string> header;

    // Check if the file is open
    if (!file.is_open())
    {
        std::cerr << "Could not open the file!" << std::endl;
        return header;
    }

    std::string line;
    if (std::getline(file, line))
    { // Read the first line
        std::stringstream ss(line);
        std::string cell;

        // Split the line by commas
        while (std::getline(ss, cell, ','))
        {
            header.push_back(cell); // Add each value to the header vector
        }
    }

    file.close();
    return header;
}

std::vector<std::vector<std::string>> getCSVRows(const std::string &filePath)
{
    std::ifstream file(filePath);
    std::vector<std::vector<std::string>> rows;

    // Check if the file is open
    if (!file.is_open())
    {
        std::cerr << "Could not open the file!" << std::endl;
        return rows;
    }

    std::string line;

    // Skip the header line
    if (std::getline(file, line))
    {
        // Continue to the next lines for rows
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            std::vector<std::string> row;

            // Split the line by commas
            while (std::getline(ss, cell, ','))
            {
                row.push_back(cell); // Add each cell to the row vector
            }
            if (row.size() > 1)
                rows.push_back(row); // Add the row to the rows vector
        }
    }

    file.close();
    return rows;
}

int main()
{
    const string ref_gnome_file_path = "/home/tarunpal/Desktop/temp/Mayank_Project/EPI_ISL_402124-ORF1ab.fasta",
                 processed_data_csv_file_path = "/home/tarunpal/Desktop/temp/Mayank_Project/processed_bam_data.csv";
    string ref_gnome = getSecondLine(ref_gnome_file_path);
    int len_ref_genome = ref_gnome.size();
    // vector<string> header = getCSVHeader(processed_data_csv_file_path);
    // vector<vector<string>> file_rows = getCSVRows(processed_data_csv_file_path);
    vector<vector<int>> processedReadData(len_ref_genome, vector<int>(2, 0));
    // string match = "Matched", unmatch = "Unmatched";
    // for (int i=0;i<file_rows.size();i++)
    // {
    //     auto row = file_rows[i];
    //     string &read = row[Column::Read];
    //     int start_idx = stoi(row[Column::Start_Index],nullptr,10);
    //     int len_read = read.size();
    //     // cout<<"Row number : "<<i+1<<" "<<read<<" \n";
    //     for (int i = 0; i < len_read; i++)
    //     {
    //         if (ref_gnome[start_idx + i] == read[i])
    //             processedReadData[start_idx + i][match]++;
    //         else
    //             processedReadData[start_idx + i][unmatch]++;
    //     }
    // }

    std::ifstream file(processed_data_csv_file_path);
    vector<std::string> row;
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Could not open the file!" << std::endl;
        return 0;
    }
    if (!std::getline(file, line)) {
        throw("Given Processed BAM File is Empty...!");
    }
    while (std::getline(file, line))
    {
        row.clear();
        std::stringstream ss(line);
        std::string cell;

        // Split the line into columns
        while (std::getline(ss, cell, ','))
        {
            row.push_back(cell);
        }

        // Process the row if it has sufficient columns
        if (row.size() > 1)
        {
            // cout<<row[Column::Read]<<"\n";
            int start_idx = std::stoi(row[Column::Start_Index],nullptr,10);
            const std::string &read = row[Column::Read];

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
                        processedReadData[position][0]++; // match++
                    }
                    else
                    {
                        processedReadData[position][1]++; // unmatch++
                    }
                }
            }
        }
    }

    file.close();

    cout << "Data Processed Sucessfully\n";

    return 0;
}