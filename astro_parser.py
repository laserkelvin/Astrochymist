
""" astro_parse.py

Script that will parse the Astrochymist bibliography of molecules
and dump the year, source, molecule, and observation technique to
a CSV file.

Current issues:
    There are no easy ways to separate the molecular formula
    and name from the way the HTML is written. If David
    ever gets around to fixing the HTML, then we can talk about
    getting this done properly.
"""

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
import requests


def main():
    url = "http://www.astrochymist.org/astrochymist_ism.html"
    response = requests.get(url)

    soup = BeautifulSoup(response.text, 'lxml')  # Parse the HTML as a string

    table = soup.find_all('table')[2]            # Grab the third table, which contains the molecules

    data = list()

    colormapping = {
        "cyan": "Radio",
        "yellow": "UV/Vis",
        "pink": "IR"
    }

    for row_index, row in enumerate(table.find_all("tr")):
        # Skip the first row
        if row_index > 0:
            # This condition is unnecessary; if you are testing
            # then drop the maximum row
            if row_index < 1000:
                columns = row.find_all(["th", "td"])
                clean_row = list()
                # Cleaning and preprocessing
                row_text = [td.text.split("\n") for td in columns]
                for col in row_text:
                    clean_row.append([x.replace("  ", "") for x in col if x])
                # Skip rows if we find data that shouldn't be in that row
                if clean_row[1][0].isdigit() or clean_row[1][0] == '4-13-1' or len(clean_row) < 4:
                    pass
                else:
                    # row_data contains all of the parsed stuff as a list
                    row_data = list()
                    for col_index, col in enumerate(clean_row):
                        if col_index == 0:
                            # Check that the first element is a number
                            # to ensure it's the year we're parsing
                            if col[0].isdigit():
                                row_data.append(int(col[0]))
                        # Second column are the molecules
                        elif col_index == 1:
                            row_data.append(col[0])
                        # Fourth column contains the source and technique
                        elif col_index == 3:
                            techs = list()
                            # Reference the last column
                            last_col = columns[-1]
                            row_data.append(col[0])
                            # Check for the color of the text, which reflects
                            # the method used for detection
                            for color, tech in colormapping.items():
                                if last_col.find("font", {"color": color}):
                                    techs.append(tech)
                            row_data.append(techs)
                    data.append(row_data)

    # Set up pandas dataframe
    df = pd.DataFrame(data=data, columns=["Year", "Molecule", "Source", "Detection Method"])

    # Remove trans-HCOOCH3 until David figures out what he wants to do
    df = df[df["Molecule"] != "trans HCOOCH3"]

    # Dump to file
    df.to_csv("astrochymist.csv", index=False)
    return df
    
if __name__ == "__main__":
    df = main()
