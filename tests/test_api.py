import inspect

from finder_charts import Crosshair, Survey, Target, create_finder_charts


def test_public_api_is_exported():
    assert Crosshair.CIRCLE_DOT.value == 1
    assert Survey.GAIA.value == "Gaia"
    assert callable(create_finder_charts)


def test_target_preserves_supplied_metadata():
    target = Target("Gaia DR3", 2016.0, 12.3, -45.6, marker_size=7, marker_color="blue", survey=Survey.GAIA)

    assert (target.catalog, target.epoch, target.ra, target.dec) == ("Gaia DR3", 2016.0, 12.3, -45.6)
    assert (target.marker_size, target.marker_color, target.survey) == (7, "blue", Survey.GAIA)


def test_chart_function_retains_its_required_coordinates():
    signature = inspect.signature(create_finder_charts)

    assert list(signature.parameters)[:2] == ["ra", "dec"]
    assert signature.parameters["img_size"].default == 100
    assert signature.parameters["legacy_surveys_timeseries"].default is False
