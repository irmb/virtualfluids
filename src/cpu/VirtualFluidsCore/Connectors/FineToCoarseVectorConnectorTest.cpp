#include <gmock/gmock.h>

#include "FineToCoarseVectorConnector.h"
#include "Block3D.h"
#include "CbVector.h"
#include "CreateTransmittersHelper.h"

class FineToCoarseVectorConnectorTest : public testing::Test
{

    void SetUp() override
    {
        block = std::make_shared<Block3D>();
    }

    void TearDown() override
    {
    
    }

public:
    CreateTransmittersHelper::TransmitterPtr senderFCevenEvenSW, receiverFCevenEvenSW;
    std::shared_ptr<Block3D> block;
};

TEST_F(FineToCoarseVectorConnectorTest, getLocalMinMax)
{
    int sendDir = D3Q27System::DIR_P00;
    block->setInterpolationFlagFC(sendDir);
    //FineToCoarseVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
                                //int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);
    InterpolationProcessorPtr iprocessor;
    auto sut = FineToCoarseVectorConnector<TbTransmitter<CbVector<LBMReal>>>(block, senderFCevenEvenSW, receiverFCevenEvenSW, sendDir, iprocessor, EvenOddNW);


    //(int &minX1, int &minX2, int &minX3, int &maxX1, int &maxX2, int &maxX3);
    //SPtr<DistributionArray3D> fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
    int maxX1 = 5; //(int)fFrom->getNX1();
    int maxX2 = 5; //(int)fFrom->getNX2();
    int maxX3 = 5;//(int)fFrom->getNX3();
    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;
    sut.getLocalMinMax(minX1, minX2, minX3, maxX1, maxX2, maxX3);

    int expectedMaxX1 = 2;
    EXPECT_THAT(maxX1, testing::Eq(expectedMaxX1));
}