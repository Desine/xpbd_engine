
load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

# --------------------------------------------------
# SFML 3.1.0
# --------------------------------------------------

cc_library(
    name = "sfml",
    linkopts = [
        "-lsfml-graphics",
        "-lsfml-window",
        "-lsfml-system",
        "-lsfml-audio",
        "-lsfml-network",
        "-lGL",
        "-fopenmp",
        "-pthread",
        "-ldl",
    ],
    visibility = ["//visibility:public"],
)

# --------------------------------------------------
# Project headers
# --------------------------------------------------

cc_library(
    name = "project_headers",
    hdrs = glob(
        [
            "include/**/*.h",
            "include/**/*.hpp",
            "src/*.h",
            "src/*.hpp",
        ],
        allow_empty = True,
    ),
    includes = [
        "include",
        "src",
    ],
    visibility = ["//visibility:public"],
)

# --------------------------------------------------
# XPBD Engine
# --------------------------------------------------

cc_binary(
    name = "exe",
    srcs = glob([
        "src/*.cpp",
    ]),
    deps = [
        ":sfml",
        ":project_headers",

        "@glm_src//:glm",
        "@imgui_src//:imgui",
        "@imgui_sfml_src//:imgui_sfml",
        "@tracy_src//:tracy",
        "@nlohmann_json//:json",
    ],
    copts = [
        "-std=c++26",
        "-fopenmp",
    ],
)
