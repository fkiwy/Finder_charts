from finder_charts import wsa


def test_extract_url_returns_download_link():
    html = '<html><a class="dl" href="https://example.test/result.csv">download</a></html>'

    assert wsa.extract_url(html) == "https://example.test/result.csv"


def test_extract_url_returns_none_without_download_link():
    assert wsa.extract_url("<html><body>No result</body></html>") is None
