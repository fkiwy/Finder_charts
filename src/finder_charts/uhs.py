import requests
from astropy.table import Table
from bs4 import BeautifulSoup


def download_html(ra, dec, radius, database):
    """
    Downloads HTML content from the WSA database for a given region.

    Parameters:
        ra (str): Right Ascension coordinate in degrees.
        dec (str): Declination coordinate in degrees.
        radius (int): Radius of the region in arcminutes.
        database (str): Name of the database (e.g., 'UHSDR3').

    Returns:
        str: HTML content if the download is successful, None otherwise.
    """
    base_url = "http://wsa.roe.ac.uk:8080/wsa/WSASQL"
    params = {
        "database": database,
        "programmeID": 107,
        "from": "source",
        "formaction": "region",
        "ra": ra,
        "dec": dec,
        "sys": "J",
        "radius": radius / 60,
        "xSize": "",
        "ySize": "",
        "boxAlignment": "RADec",
        "format": "CSV",
        "select": "default",
        "where": "",
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to download HTML. Status code: {response.status_code}")
        print(response.text)
        return None


def extract_url(html_content):
    """
    Extracts the CSV file URL from the HTML content.

    Parameters:
        html_content (str): HTML content obtained from the WSA database.

    Returns:
        str: File URL if found, None otherwise.
    """
    soup = BeautifulSoup(html_content, "html.parser")
    download_link = soup.find("a", class_="dl")
    if download_link:
        return download_link["href"]
    else:
        return None


def download_file(url):
    """
    Downloads and converts the CSV file to an Astropy table.

    Parameters:
        url (str): URL of the CSV file.

    Returns:
        astropy.table.Table: Astropy table if successful, None otherwise.
    """
    response = requests.get(url)
    if response.status_code == 200:
        try:
            return Table.read(
                response.text,
                format="ascii",
                comment="#"
            )
        except IndexError:
            return None
    else:
        print(f"Failed to download CSV file. Status code: {response.status_code}")
        print(response.text)
        return None


def query_region(ra, dec, radius=5, database="UHSDR3"):
    """
    Queries the WSA database for a specific region and returns the results as an Astropy table.

    Parameters:
        ra (str): Right Ascension coordinate in degrees.
        dec (str): Declination coordinate in degrees.
        radius (int): Radius of the region in arcminutes (default is 5 arcminutes).
        database (str): Name of the database (default is 'UHSDR3').

    Returns:
        astropy.table.Table: Astropy table containing the query results, or None if unsuccessful.
    """
    table = None
    html_content = download_html(ra, dec, radius, database)
    if html_content:
        url = extract_url(html_content)
        if url:
            table = download_file(url)
            if table:
                table.sort("distance")
    return table
