"""Fail closed if a committed scientific input changes."""

from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "tess2019357164649-s0020-0000000086396382-0165-s_lc.fits": "497b2e27ffc8f6a2f1c76b388d7a5f4943ecd4e48b89c3882e7d014be91857b7",
    "tess2021258175143-s0043-0000000086396382-0214-s_lc.fits": "e3d5c6f76884014ec1dec1eb60f04eaf734330a7abd8eef99141f49df03b2493",
    "tess2021284114741-s0044-0000000086396382-0215-s_lc.fits": "fab2038adf1134d3d5069c2dc5e01f2bcd17a7d524629edb669118a11b204e0a",
    "tess2021310001228-s0045-0000000086396382-0216-s_lc.fits": "1ac919f1a4b472b7d6c1911901e4f967e222cec8559204c58c3b091d995f5d06",
    "tess2023289093419-s0071-0000000086396382-0266-s_lc.fits": "2e112a9352bab08e4e3228daebb9b53d6fe829378fea7b5a4aafcbed7187e844",
    "tess2023315124025-s0072-0000000086396382-0267-s_lc.fits": "5330ed8251ef5ece4c96ecb55023652f1014a2a007af4df9fcc438fd9fbdb586",
    "published_transit_occultation_times.csv": "256924da97037ed037b3cf327ee7206d9099382d6ab215dff134c09c4351f3a1",
    "leonardi2024_literature_timings.tsv": "14688a4142526a54d60a9f872f159506f5584befb9d1fa4ce6f6366b24fc2e39",
}


def test_scientific_inputs_match_documented_sha256():
    actual_files = {
        path.name for path in (ROOT / "data").glob("tess*_lc.fits")
    } | {"published_transit_occultation_times.csv", "leonardi2024_literature_timings.tsv"}
    assert actual_files == set(EXPECTED)
    for filename, expected in EXPECTED.items():
        path = ROOT / "data" / filename
        payload = (
            path.read_text(encoding="utf-8").replace("\r\n", "\n").encode()
            if path.suffix == ".csv"
            else path.read_bytes()
        )
        digest = sha256(payload).hexdigest()
        assert digest == expected, f"checksum mismatch: {filename}"
