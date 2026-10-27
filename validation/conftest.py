def pytest_addoption(parser):
    parser.addoption(
        "--opengate",
        choices=["json", "images"],
        default="json",
        help="OpenGATE validation level: 'json' compares with the stored OpenGATE values "
        "(default, no download); 'images' uses the real images and OpenGATE outputs "
        "(downloaded from Zenodo on first use).",
    )
