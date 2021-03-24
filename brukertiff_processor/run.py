import os

import bioformats
import javabridge

from brukertiff_processor import BRUKERTIFFProcessor

if __name__ == '__main__':
    try:
        log_config = os.path.join(os.path.split(__file__)[0], "resources", "log4j.properties")

        javabridge.start_vm(
            args=[
                "-Dlog4j.configuration=file:{}".format(log_config),
            ],
            class_path=bioformats.JARS,
            run_headless=True,
            max_heap_size="4096m"
        )

        myloglevel = "ERROR"  # user string argument for logLevel.
        rootLoggerName = javabridge.get_static_field("org/slf4j/Logger", "ROOT_LOGGER_NAME", "Ljava/lang/String;")
        rootLogger = javabridge.static_call("org/slf4j/LoggerFactory", "getLogger",
                                            "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
        logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level", myloglevel,
                                               "Lch/qos/logback/classic/Level;")
        #javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)

        task = BRUKERTIFFProcessor(cli=True)
        task.run()
    finally:
        javabridge.kill_vm()
