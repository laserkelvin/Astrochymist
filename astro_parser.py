
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
import re


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
            if row_index < 400:
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
                        elif col_index == 2:
                            # Loop over the references, and use regex to extract
                            # the author names
                            authors = list()
                            for text in col:
                                # Regex will pick up capital letter followed by period,
                                # space, a capital letter and at least three lower case
                                authors.extend(
                                    re.findall(r"[A-Z]\.\s[A-Za-z]\w{3,}", text)
                                )
                                # Find all authors not formatted the same way
                                authors.extend(
                                    re.findall(r"[A-Z]\.[A-Za-z]\w{3,}", text)
                                )
                            # Remove whitespace
                            authors = [author.replace(" ", "") for author in authors]
                            # Make sure all authors are formatted the same way
                            authors = [" ".join(author.split(".")) for author in authors]
                            row_data.append(authors)
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
    df = pd.DataFrame(data=data, columns=["Year", "Molecule", "Authors", "Source", "Detection Method"])

    # Remove trans-HCOOCH3 until David figures out what he wants to do
    df = df[df["Molecule"] != "trans HCOOCH3"]

    # Dump to file
    df.to_csv("astrochymist.csv", index=False)
    return df

def time_analysis(df):
    # Perform analysis of detections of time, as well as
    # the number of discoveries per author over time.
    author_df = pd.DataFrame(df["Authors"].values.tolist(), index=df["Year"])
    unique_authors = list()
    # Create a list of unique authors
    for column in author_df:
        for author in author_df[column].unique():
            if author not in unique_authors:
                unique_authors.append(author)
    # Dictionary for holding stuff
    authorship = list()
    tally = list()
    running_tally = {author: 0 for author in unique_authors}
    for index, row in df.iterrows():
        row_data = {author: False for author in unique_authors}
        for author in row["Authors"]:
            row_data[author] = True
        authorship.append(row_data)
    authorship_df = pd.DataFrame(authorship, index=df["Year"])
    # Collapse the dataframe along years
    authorship_df = authorship_df.groupby(authorship_df.index).sum()
    # Cumulative sum of authors over the years
    tally_df = authorship_df.cumsum(axis=0)
    tally_df.index = authorship_df.index
    # Get the ranking of detections by authors
    rank_df = authorship_df.sum().sort_values(ascending=False) 
    authorship_df.to_csv("author_timeseries.csv")
    tally_df.to_csv("author_cumsum.csv")
    rank_df.to_csv("author_ranking.csv")
    return authorship_df, tally_df, rank_df
    
if __name__ == "__main__":
    df = main()
