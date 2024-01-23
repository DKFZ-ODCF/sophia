#ifndef FIXTURES_H_
#define FIXTURES_H_

#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <optional>
#include <memory>
#include "global.h"
#include "ChrInfoTable.h"
#include "Hg38ChrConverter.h"
#include "GlobalAppConfig.h"


namespace sophia {

    class ChrInfoTableFixture : public ::testing::Test {

      protected:

        ChrInfoTable *chr_info_table = nullptr;

        void SetUp() override {
            std::string chromosome_file = "resources/hg38_test.tsv";
            std::vector<ChrInfo> chr_info = read_chr_info(chromosome_file);
            chr_info_table = new ChrInfoTable(chr_info);
        }

        void TearDown() override {
            delete chr_info_table;
        }

    };

    class Hg38ChrConverterFixture : public ChrInfoTableFixture {

      protected:

        void SetUp() {
            ChrInfoTableFixture::SetUp();
            try {
                std::unique_ptr<Hg38ChrConverter> converter =
                    std::unique_ptr<Hg38ChrConverter>(new Hg38ChrConverter("hg38", *chr_info_table));
                // TODO Find a way to deal with the singleton, without having to refactor the whole code
                GlobalAppConfig::init(std::move(converter));
            } catch (const std::logic_error& e) {
                // In case the singleton is already set, this will throw a logic_error, which we
                // we just ignore.
            }
        }

        void TearDown() override {
            // TODO Maybe remove the singleton here? At least for non-parallel testing
        }

    };

}

#endif /* FIXTURES_H_ */
