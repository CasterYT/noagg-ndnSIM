//
// Created by 李怿曈 on 16/4/2024.
//

#include "ndn-aggregator.hpp"
#include "model/ndn-l3-protocol.hpp"
#include "helper/ndn-fib-helper.hpp"
#include "ModelData.hpp"
#include "ndn-consumer.hpp"

#include "ns3/ptr.h"
#include "ns3/log.h"
#include "ns3/simulator.h"
#include "ns3/packet.h"
#include "ns3/callback.h"
#include "ns3/string.h"
#include "ns3/boolean.h"
#include "ns3/uinteger.h"
#include "ns3/integer.h"
#include "ns3/double.h"
#include <limits>

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>
#include <iostream>
#include <sstream>

#include "utils/ndn-ns3-packet-tag.hpp"
#include "utils/ndn-rtt-mean-deviation.hpp"

#include <ndn-cxx/lp/tags.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/ref.hpp>

NS_LOG_COMPONENT_DEFINE("ndn.Aggregator");

namespace ns3{
namespace ndn{

NS_OBJECT_ENSURE_REGISTERED(Aggregator);

TypeId
Aggregator::GetTypeId(void)
{
    static TypeId tid =
            TypeId("ns3::ndn::Aggregator")
            .SetGroupName("Ndn")
            .SetParent<App>()
            .AddConstructor<Aggregator>()
            .AddAttribute("StartSeq", "Starting sequence number", IntegerValue(0),
                          MakeIntegerAccessor(&Aggregator::m_seq), MakeIntegerChecker<int32_t>())
            .AddAttribute("Prefix", "Interest prefix/name", StringValue("/"),
                          MakeNameAccessor(&Aggregator::m_prefix), MakeNameChecker())
            .AddAttribute("NodePrefix", "Node prefix", StringValue(),
                          MakeStringAccessor(&Aggregator::m_nodeprefix), MakeStringChecker())
            .AddAttribute("LifeTime", "Life time for interest packet", StringValue("4s"),
                          MakeTimeAccessor(&Aggregator::m_interestLifeTime), MakeTimeChecker())
            .AddAttribute("RetxTimer",
                          "Timeout defining how frequent retransmission timeouts should be checked",
                          StringValue("50ms"),
                          MakeTimeAccessor(&Aggregator::GetRetxTimer, &Aggregator::SetRetxTimer),
                          MakeTimeChecker())
            .AddAttribute("Freshness", "Freshness of data packets, if 0, then unlimited freshness",
                          TimeValue(Seconds(0)), MakeTimeAccessor(&Aggregator::m_freshness),
                          MakeTimeChecker())
            .AddAttribute("Signature","Fake signature, 0 valid signature (default), other values application-specific",
                          UintegerValue(0), MakeUintegerAccessor(&Aggregator::m_signature),
                          MakeUintegerChecker<uint32_t>())
            .AddAttribute("KeyLocator",
                          "Name to be used for key locator.  If root, then key locator is not used",
                          NameValue(), MakeNameAccessor(&Aggregator::m_keyLocator), MakeNameChecker())
            .AddAttribute("Window", "Initial size of the window", StringValue("1"),
                          MakeUintegerAccessor(&Aggregator::GetWindow, &Aggregator::SetWindow),
                          MakeUintegerChecker<uint32_t>())
            .AddAttribute("MaxSeq", "Maximum sequence number to request (alternative to Size attribute, "
                                    "would activate only if Size is -1). "
                                    "The parameter is activated only if Size negative (not set)",
                            IntegerValue(std::numeric_limits<uint32_t>::max()),
                            MakeUintegerAccessor(&Aggregator::GetSeqMax, &Aggregator::SetSeqMax),
                            MakeUintegerChecker<uint32_t>())
            .AddAttribute("InitialWindowOnTimeout", "Set window to initial value when timeout occurs",
                          BooleanValue(true),
                          MakeBooleanAccessor(&Aggregator::m_setInitialWindowOnTimeout),
                          MakeBooleanChecker())
            .AddAttribute("Beta", "TCP Multiplicative Decrease factor", DoubleValue(0.5),
                          MakeDoubleAccessor(&Aggregator::m_beta),
                          MakeDoubleChecker<double>())
            .AddAttribute("AddRttSuppress", "Minimum number of RTTs (1 + this factor) between window decreases",
                          DoubleValue(0.5),
                          MakeDoubleAccessor(&Aggregator::m_addRttSuppress),
                          MakeDoubleChecker<double>())
            .AddAttribute("ReactToCongestionMarks",
                          "If true, process received congestion marks",
                          BooleanValue(true),
                          MakeBooleanAccessor(&Aggregator::m_reactToCongestionMarks),
                          MakeBooleanChecker())
            .AddAttribute("UseCwa",
                          "If true, use Conservative Window Adaptation",
                          BooleanValue(true),
                          MakeBooleanAccessor(&Aggregator::m_useCwa),
                          MakeBooleanChecker())
            .AddTraceSource("LastRetransmittedInterestDataDelay",
                            "Delay between last retransmitted Interest and received Data",
                            MakeTraceSourceAccessor(&Aggregator::m_lastRetransmittedInterestDataDelay),
                            "ns3::ndn::Aggregator::LastRetransmittedInterestDataDelayCallback")
            .AddTraceSource("FirstInterestDataDelay",
                            "Delay between first transmitted Interest and received Data",
                            MakeTraceSourceAccessor(&Aggregator::m_firstInterestDataDelay),
                            "ns3::ndn::Aggregator::FirstInterestDataDelayCallback")
            .AddTraceSource("WindowTrace",
                            "Window that controls how many outstanding interests are allowed",
                            MakeTraceSourceAccessor(&Aggregator::m_window),
                            "ns3::ndn::Aggregator::WindowTraceCallback")
            .AddTraceSource("InFlight", "Current number of outstanding interests",
                            MakeTraceSourceAccessor(&Aggregator::m_inFlight),
                            "ns3::ndn::Aggregator::WindowTraceCallback");

    return tid;

}

Aggregator::Aggregator()
    : m_rand(CreateObject<UniformRandomVariable>())
    , m_inFlight(0)
    , m_ssthresh(std::numeric_limits<double>::max())
    , m_highData(0)
    , m_recPoint(0.0)
    , m_seq(0)
    , totalAggregateTime(0)
    , totalResponseTime(0)
    , round(0)
    , iteration(0)
{

    m_rtt = CreateObject<RttMeanDeviation>();
    //m_timeoutThreshold = ns3::Seconds(0.3);
}




// Helper function to split string by a delimiter and return vector of strings
std::pair<std::string, std::set<std::string>>
Aggregator::aggTreeProcessSingleString(const std::string& input)
{
    std::istringstream iss(input);
    std::string segment;
    std::vector<std::string> segments;

    // Use getline to split the string by '.'
    while (getline(iss, segment, '.')) {
        segments.push_back(segment);
    }

    // Check if there are enough segments to form a key and a set
    if (segments.size() > 1) {
        std::string key = segments[0];
        std::set<std::string> values(segments.begin() + 1, segments.end());
        return {key, values};
    }

    return {};  // Return an empty pair if not enough segments
}


std::map<std::string, std::set<std::string>>
Aggregator::aggTreeProcessStrings(const std::vector<std::string>& inputs)
{
    std::map<std::string, std::set<std::string>> result;

    for (const std::string& input : inputs) {
        auto entry = aggTreeProcessSingleString(input);
        if (!entry.first.empty()) {
            result[entry.first].insert(entry.second.begin(), entry.second.end());
        }
    }

    return result;
}




std::map<std::string, std::set<std::string>>
Aggregator::getLeafNodes(const std::string& key,
            const std::map<std::string, std::vector<std::string>>& treeMap)
{
    return App::getLeafNodes(key, treeMap);
}


void
Aggregator::ResponseTimeSum (int64_t response_time)
{
    totalResponseTime += response_time;
    ++round;
}

int64_t
Aggregator::GetResponseTimeAverage()
{
    if (round == 0)
    {
        NS_LOG_DEBUG("Error happened when calculating average response time!");
        return 0;
    }

    return totalResponseTime / round;
}


void
Aggregator::AggregateTimeSum (int64_t aggregate_time)
{
    totalAggregateTime += aggregate_time;
    ++iteration;
}

int64_t
Aggregator::GetAggregateTimeAverage()
{
    if (iteration == 0)
    {
        NS_LOG_DEBUG("Error happened when calculating aggregate time!");
        return 0;
    }

    return totalAggregateTime / iteration;
}



// New design to get child info's mapping
void
Aggregator::AddLinkInfo(const std::string& filename)
{
    App::AddLinkInfo(filename);

    std::ifstream file(filename);
    std::string line;
    bool parseSection = false;

    if (!file.is_open()) {
        NS_LOG_ERROR("Unable to open file " << filename);
        return;
    }

    while (getline(file, line)) {
        // Check for section header
        if (line.find("link") != std::string::npos) {
            parseSection = true;
            continue;
        }

        if (parseSection && !line.empty()) {
            std::istringstream iss(line);
            std::string parentNode, childNode;

            // Read the first two columns
            if (iss >> parentNode >> childNode) {
                // store parent-child mapping into global variable
                m_linkInfo[parentNode].push_back(childNode);
            }
        }
    }

    file.close();
}


void
Aggregator::CheckRetxTimeout()
{
    Time now = Simulator::Now();

    Time rto = m_rtt->RetransmitTimeout();

    while (!m_seqTimeouts.empty()) {
        SeqTimeoutsContainer::index<i_timestamp>::type::iterator entry =
                m_seqTimeouts.get<i_timestamp>().begin();
        if (entry->time + rto <= now)
        {
            uint32_t seqNo = entry->seq;
            m_seqTimeouts.get<i_timestamp>().erase(entry);
            OnTimeout(seqNo);
        } else
            break;
    }

    m_retxEvent = Simulator::Schedule(m_retxTimer, &Aggregator::CheckRetxTimeout, this);

/*    for (auto it = m_timeoutCheck.begin(); it != m_timeoutCheck.end();){
        if (now - it->second > m_timeoutThreshold) {
            uint32_t seq = it->first;
            OnTimeout(seq);
            it = m_timeoutCheck.erase(it);
        } else {
            ++it;
        }
    }
    Simulator::Schedule(Seconds(0.05), &Aggregator::CheckRetxTimeout, this);*/
}

void
Aggregator::OnTimeout(uint32_t sequenceNumber)
{
    NS_LOG_FUNCTION(sequenceNumber);

    m_rtt->IncreaseMultiplier(); // Double the next RTO
    m_rtt->SentSeq(SequenceNumber32(sequenceNumber),
                   1); // make sure to disable RTT calculation for this sample
    m_retxSeqs.insert(sequenceNumber);

    /// Designed for AIMD
    WindowDecrease();

    if (m_inFlight > static_cast<uint32_t>(0)){
        m_inFlight--;
    }
    NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);


    // Customized design, when timeout happens, retrieve vector (consisting the names related to this sequence number) and names stored inside for retransmission
    if (m_newSeq_Name.find(sequenceNumber) != m_newSeq_Name.end()) {
        for (auto item : m_newSeq_Name[sequenceNumber]) {
            shared_ptr<Name> name = make_shared<Name>(item);
            NS_LOG_INFO("Retransmitting interest: " << sequenceNumber);
            SendInterest(name, sequenceNumber);
        }
    } else {
        NS_LOG_INFO("No mapping found for sequence number: " << sequenceNumber);
    }

}


void
Aggregator::SetRetxTimer(Time retxTimer)
{
    m_retxTimer = retxTimer;
    if (m_retxEvent.IsRunning()) {
        Simulator::Remove(m_retxEvent);
    }

    // Schedule new timeout
    m_retxEvent = Simulator::Schedule(m_retxTimer, &Aggregator::CheckRetxTimeout, this);
}

Time
Aggregator::GetRetxTimer() const
{
    return m_retxTimer;
}

void
Aggregator::StartApplication()
{
    NS_LOG_FUNCTION_NOARGS();
    App::StartApplication();
    FibHelper::AddRoute(GetNode(), m_prefix, m_face, 0);

    //AddLinkInfo(fileAggBottleneck);
    //Simulator::Schedule(Seconds(0.05), &Aggregator::CheckRetxTimeout, this);
}

void
Aggregator::StopApplication()
{
    NS_LOG_FUNCTION_NOARGS();
    // Cancel packet generation
    Simulator::Cancel(m_sendEvent);

    NS_LOG_INFO("The average response time of Aggregator in " << round << " aggregation rounds is: " << GetResponseTimeAverage() << " ms");
    NS_LOG_INFO("The average aggregate time is: " << GetAggregateTimeAverage() << " ms");
    App::StopApplication();
}

void Aggregator::aggregate(const ModelData& data, const std::string& dataName) {
    // first initialization
    if (sumParameters.find(dataName) == sumParameters.end()){
        sumParameters[dataName] = std::vector<float>(300, 0.0f);
        count[dataName] = 0;
    }

    // Aggregate data
    std::transform(sumParameters[dataName].begin(), sumParameters[dataName].end(), data.parameters.begin(), sumParameters[dataName].begin(), std::plus<float>());
    count[dataName]++;

/*    // Print the data from data.parameters
    std::cout << "Aggregator's aggregate function processed Parameters for '" << dataName << "':" << std::endl;
    for (size_t i = 0; i < data.parameters.size(); ++i) {
        std::cout << "Parameter " << i << ": " << data.parameters[i] << std::endl;
    }*/
}

ModelData Aggregator::getMean(const std::string& dataName){
    ModelData result;
    if (sumParameters.find(dataName) != sumParameters.end()) {
        result.parameters = sumParameters[dataName];  // Direct assignment
    }

/*    // Print the data from data.parameters
    std::cout << "Aggregator getMean function processed Parameters for '" << dataName << "':" << std::endl;
    for (size_t i = 0; i < result.parameters.size(); ++i) {
        std::cout << "Parameter " << i << ": " << result.parameters[i] << std::endl;
    }*/

    return result;
}

void
Aggregator::WillSendOutInterest(uint32_t sequenceNumber)
{
    //NS_LOG_DEBUG("Trying to add " << sequenceNumber << " with " << Simulator::Now() << ". already " << m_seqTimeouts.size() << " items");

    // Start calculating Timeouts
    m_seqTimeouts.insert(SeqTimeout(sequenceNumber, Simulator::Now()));
    // Delay for certain sequence number
    m_seqFullDelay.insert(SeqTimeout(sequenceNumber, Simulator::Now()));

    m_seqLastDelay.erase(sequenceNumber);
    m_seqLastDelay.insert(SeqTimeout(sequenceNumber, Simulator::Now()));
    // Retransmission count
    m_seqRetxCounts[sequenceNumber]++;
    m_rtt->SentSeq(SequenceNumber32(sequenceNumber), 1);

    //m_timeoutCheck[sequenceNumber] = Simulator::Now();

    /// Designed for Window
    m_inFlight++;
}

void
Aggregator::OnNack(shared_ptr<const lp::Nack> nack)
{
    /// tracing inside
    App::OnNack(nack);

    NS_LOG_INFO("NACK received for: " << nack->getInterest().getName()
                                          << ", reason: " << nack->getReason());
}

void
Aggregator::OnInterest(shared_ptr<const Interest> interest)
{
    NS_LOG_INFO("Receiving interest:  " << *interest);
    NS_LOG_INFO("The incoming interest packet size is: " << interest->wireEncode().size());
    App::OnInterest(interest);

    std::string interestType = interest->getName().get(-2).toUri();



    NS_LOG_INFO("interestType: " << interestType);
    if (interestType == "data") {

        // Parse incoming interest, retrieve their name segments, currently use "/NextHop/Destination/Type/SeqNum"
        //std::string parent_prefix = m_prefix.toUri().substr(1);
        std::string nexthop = interest->getName().get(0).toUri();
        std::string dest = interest->getName().get(1).toUri();
        std::string seqNum = interest->getName().get(-1).toUri();
        uint32_t seq = interest->getName().get(-1).toSequenceNumber();

        std::vector<std::string> segments;  // Vector to store the segments
        std::istringstream iss(dest);
        std::string segment;
        std::vector<std::string> seqName; // vector to store name, need to be stored in the mapping, used for retransmission later

        // design for data aggregation, store the info about how interest is divided into segments
        std::vector<std::string> value_agg;

        // split the destination segment into several ones and store them individually in a vector called "segments"
        while (std::getline(iss, segment, '.')) {
            segments.push_back(segment);
        }

        /// get current time for aggregation time calculation, need to change later!!!!!!!!!!!!!
        aggregateStartTime[seq] = ns3::Simulator::Now();


        // iterate all child nodes in a loop, assign the "name" based on child node's prefix
        for (const auto& [child, leaves] : aggregationMap) {
            std::string name_sec1;
            std::string name;

            // interest is divided
            for (const auto& leaf : leaves) {
                if (std::find(segments.begin(), segments.end(), leaf) != segments.end()) {
                    name_sec1 += leaf + ".";
                } else {
                    NS_LOG_INFO("Data from " << leaf << " is not required for this iteration.");
                }
            }

            if (name_sec1.empty()) {
                NS_LOG_INFO("No interest needs to be sent to node " << child << " in this iteration");
                continue;
            } else {
                name_sec1.resize(name_sec1.size() - 1);
                value_agg.push_back(name_sec1);

                /// for calculating response time, need to change later
                std::string resKey = "/" + name_sec1 + "/" + seqNum;
                currentTime[resKey] = ns3::Simulator::Now();

                name = "/" + child + "/" + name_sec1 + "/data";
                shared_ptr<Name> newName = make_shared<Name>(name);
                seqName.push_back(name);

                /// Store relevant into FIFO queue, send later!!!!!!!!!!!
                interestQueue.push(std::make_tuple(newName, seq, resKey));
            }

        }

        if (map_agg_oldSeq_newName.find(seqNum) == map_agg_oldSeq_newName.end()){
            map_agg_oldSeq_newName[seqNum] = value_agg;
            m_agg_newDataName[seq] = dest;
/*
            // testing
            NS_LOG_DEBUG("map_agg_oldSeq_newName[" << seqNum << "]: ");
            for (const auto& item : value_agg){
                NS_LOG_DEBUG(item);
            }

            NS_LOG_DEBUG("m_agg_newDataName[" << seqNum << "]: " << dest);*/
        }

        if (m_newSeq_Name.find(seq) == m_newSeq_Name.end()){
            m_newSeq_Name[seq] = seqName; // This is wrong (seqName)!

/*            // testing
            NS_LOG_DEBUG("m_newSeq_Name[" << seqNum << "]: ");
            for (const auto& item : seqName){
                NS_LOG_DEBUG(item);
            }*/
        }

        ScheduleNextPacket();

    } else if (interestType == "initialization") {

        // Extract useful info and parse it into readable format
        std::vector<std::string> inputs;
        if (interest->getName().size() > 3) {
            for (size_t i = 1; i < interest->getName().size() - 2; ++i) {
                inputs.push_back(interest->getName().get(i).toUri());
            }
        }
        aggregationMap = aggTreeProcessStrings(inputs);


        // testing, delete later!!!!
        if (aggregationMap.empty())
            NS_LOG_DEBUG("aggregationMap is empty!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        else {
            // Print the result mapping
            for (const auto& [node, leaves] : aggregationMap) {
                NS_LOG_DEBUG(node << " contains leaf nodes: ");
                for (const auto& leaf : leaves) {
                    NS_LOG_DEBUG(leaf << " ");
                }
            }
        }

        Name dataName(interest->getName());
        auto data = make_shared<Data>();
        data->setName(dataName);
        data->setFreshnessPeriod(::ndn::time::milliseconds(m_freshness.GetMilliSeconds()));

        SignatureInfo signatureInfo(static_cast< ::ndn::tlv::SignatureTypeValue>(255));
        if (m_keyLocator.size() > 0) {
            signatureInfo.setKeyLocator(m_keyLocator);
        }
        data->setSignatureInfo(signatureInfo);
        ::ndn::EncodingEstimator estimator;
        ::ndn::EncodingBuffer encoder(estimator.appendVarNumber(m_signature), 0);
        encoder.appendVarNumber(m_signature);
        data->setSignatureValue(encoder.getBuffer());

        // to create real wire encoding
        data->wireEncode();
        m_transmittedDatas(data, this, m_face);
        m_appLink->onReceiveData(*data);

    }
}


void
Aggregator::SetWindow(uint32_t window)
{
    m_initialWindow = window;
    m_window = m_initialWindow;
}

uint32_t
Aggregator::GetWindow() const
{
    return m_initialWindow;
}

void
Aggregator::SetSeqMax(uint32_t seqMax)
{
    // Be careful, ignore maxSize here
    m_seqMax = seqMax;
}

uint32_t
Aggregator::GetSeqMax() const
{
    return m_seqMax;
}


void
Aggregator::WindowIncrease()
{
    if (m_window < m_ssthresh) {
        m_window += 1.0;
    }else {
        m_window += (1.0 / m_window);
    }
    NS_LOG_DEBUG("Window size increased to " << m_window);
}

void
Aggregator::WindowDecrease()
{
    if (!m_useCwa || m_highData > m_recPoint) {
        const double diff = m_seq - m_highData;
        BOOST_ASSERT(diff > 0);

        m_recPoint = m_seq + (m_addRttSuppress * diff);

        // AIMD
        m_ssthresh = m_window * m_beta;
        m_window = m_ssthresh;

        // Window size can't be reduced below initial size
        if (m_window < m_initialWindow) {
            m_window = m_initialWindow;
        }

        NS_LOG_DEBUG("Window size decreased to " << m_window);
    }
    else {
        NS_LOG_DEBUG("Window decrease suppressed, HighData: " << m_highData << ", RecPoint: " << m_recPoint);
    }

}


void
Aggregator::ScheduleNextPacket()
{
    if (m_window == static_cast<uint32_t>(0)) {
        Simulator::Remove(m_sendEvent);

        NS_LOG_DEBUG("New event in " << (std::min<double>(0.5, m_rtt->RetransmitTimeout().ToDouble(Time::S))) << " sec");

        m_sendEvent = Simulator::Schedule(Seconds(std::min<double>(0.5, m_rtt->RetransmitTimeout().ToDouble(Time::S))),
                                          &Aggregator::SendPacket, this);

    }
    else if (m_inFlight >= m_window) {
        // do nothing
    }
    else {
        if (m_sendEvent.IsRunning()) {
            Simulator::Remove(m_sendEvent);
        }

        m_sendEvent = Simulator::ScheduleNow(&Aggregator::SendPacket, this);
    }
}



void
Aggregator::SendPacket()
{
    if (!interestQueue.empty()){
        auto interestInfo = interestQueue.front();
        interestQueue.pop();
        shared_ptr<Name> newName = std::get<0>(interestInfo);
        uint32_t seq = std::get<1>(interestInfo);
        std::string resKey = std::get<2>(interestInfo);

        NS_LOG_INFO("Arrange packet sending: " << newName->toUri());
        SendInterest(newName, seq);

        if (currentTime.find(resKey) == currentTime.end())
            currentTime[resKey] = ns3::Simulator::Now();

        if (aggregateStartTime.find(seq) == aggregateStartTime.end())
            aggregateStartTime[seq] = ns3::Simulator::Now();

        ScheduleNextPacket();
    } else {
        NS_LOG_DEBUG("No info in the queue, aggregator can't send anything!");
    }
}

void
Aggregator::SendInterest(shared_ptr<Name> newName, uint32_t seq)
{
    if (!m_active)
        return;

    NS_LOG_INFO("Sending new interest: " << newName->toUri());
    newName->appendSequenceNumber(seq);
    shared_ptr<Interest> newInterest = make_shared<Interest>();
    newInterest->setNonce(m_rand->GetValue(0, std::numeric_limits<uint32_t>::max()));
    newInterest->setCanBePrefix(false);
    newInterest->setName(*newName);
    time::milliseconds interestLifeTime(m_interestLifeTime.GetMilliSeconds());
    newInterest->setInterestLifetime(interestLifeTime);

    WillSendOutInterest(seq);
    m_transmittedInterests(newInterest, this, m_face);
    m_appLink->onReceiveInterest(*newInterest);
}

void
Aggregator::OnData(shared_ptr<const Data> data)
{
    if(!m_active)
        return;

    App::OnData(data);
    NS_LOG_INFO ("Received content object: " << boost::cref(*data));
    NS_LOG_INFO("The incoming data packet size is: " << data->wireEncode().size());

    std::string seqNum = data->getName().get(-1).toUri();
    uint32_t seq = data->getName().at(-1).toSequenceNumber();

    // Essential operation for default hop count calculation, delay measurement, sequence number management (not needed in this case), RTT updates
    int hopCount = 0;
    auto hopCountTag = data->getTag<lp::HopCountTag>();
    if (hopCountTag != nullptr) { // e.g., packet came from local node's cache
        hopCount = *hopCountTag;
    }
    NS_LOG_DEBUG("Hop count: " << hopCount);

    SeqTimeoutsContainer::iterator entry = m_seqLastDelay.find(seq);
    if (entry != m_seqLastDelay.end()) {
        m_lastRetransmittedInterestDataDelay(this, seq, Simulator::Now() - entry->time, hopCount);
    }

    entry = m_seqFullDelay.find(seq);
    if (entry != m_seqFullDelay.end()) {
        m_firstInterestDataDelay(this, seq, Simulator::Now() - entry->time, m_seqRetxCounts[seq], hopCount);
    }

    m_seqRetxCounts.erase(seq);
    m_seqFullDelay.erase(seq);
    m_seqLastDelay.erase(seq);

    m_seqTimeouts.erase(seq);
    m_retxSeqs.erase(seq);

    m_rtt->AckSeq(SequenceNumber32(seq));

    /// Designed for AIMD
    if (m_highData < seq) {
        m_highData = seq;
    }

    if (data->getCongestionMark() > 0) {
        if (m_reactToCongestionMarks) {
            NS_LOG_DEBUG("Received congestion mark: " << data->getCongestionMark());
            WindowDecrease();
        }
        else {
            NS_LOG_DEBUG("Ignored received congestion mark: " << data->getCongestionMark());
        }
    } else {
        WindowIncrease();
    }

    if (m_inFlight > static_cast<uint32_t>(0)) {
        m_inFlight--;
    }

    NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);

    ScheduleNextPacket();


    // Check what are the exact names of data waiting for aggregation
    ModelData modelData;
    auto data_map = map_agg_oldSeq_newName.find(seqNum);
    if (data_map != map_agg_oldSeq_newName.end())
    {
        // Data name's variable definition

        NS_LOG_INFO("Received data name: " << data->getName().toUri());
        std::string data_key = "/" + data->getName().get(-1).toUri();
        //NS_LOG_INFO("data_key: " << data_key);
        std::string data_value = data->getName().get(1).toUri();
        //NS_LOG_INFO("data_value: " << data_value);
        std::string child = data->getName().get(0).toUri();

        auto& vec = data_map->second;
        auto vecIt = std::find(vec.begin(), vec.end(), data_value);
        std::vector<uint8_t> oldbuffer(data->getContent().value(), data->getContent().value() + data->getContent().value_size());

        if (deserializeModelData(oldbuffer, modelData) && vecIt != vec.end()) {
            aggregate(modelData, seqNum);
            vec.erase(vecIt);
            //map_agg_final[seqNum].push_back(data_value);
        } else{
            NS_LOG_INFO("Data name doesn't exist in map_agg_oldSeq_newName, meaning this data packet is duplicate, do nothing!");
            return;
        }

        // testing
        NS_LOG_DEBUG("Existing elements to be aggregated (map_agg_oldSeq_newName):");
        for (auto item : map_agg_oldSeq_newName[seqNum]){
            NS_LOG_DEBUG("Iterating map_agg_oldSeq_newName's item: " << item);
        }

        // for response time measurement
        std::string resKey = "/" + data_value + data_key;
        if (currentTime.find(resKey) != currentTime.end()){
            responseTime[resKey] = ns3::Simulator::Now() - currentTime[resKey];
            ResponseTimeSum(responseTime[resKey].GetMilliSeconds());
            currentTime.erase(resKey);
            //NS_LOG_INFO("Aggregator's response time of child " << child <<" and sequence " << seq << " is: " << responseTime[resKey].GetMilliSeconds() << " ms");
        }

        NS_LOG_INFO("aggregation count: " << count[seqNum]);

        // Check the mapping to judge whether the aggregation process is done
        if (vec.empty()){

            NS_LOG_DEBUG("Aggregation finished.");
            // reset the aggregation count
            count[seqNum] = 0;

            // Calculate aggregate time ！
            if (aggregateStartTime.find(seq) != aggregateStartTime.end()) {
                aggregateTime[seq] = ns3::Simulator::Now() - aggregateStartTime[seq];
                AggregateTimeSum(aggregateTime[seq].GetMilliSeconds());
                aggregateStartTime.erase(seq);
                //NS_LOG_INFO("Aggregator's aggregate time of sequence " << seq << " is: " << aggregateTime[seq].GetMilliSeconds());
            } else {
                NS_LOG_DEBUG("Error when calculating aggregation time, no reference found for seq " << seq);
            }

            // Stop checking timeout associated with this seq
            //m_timeoutCheck.erase(seq);

            // For data retransmission, erase the sequence's mapping
            auto it = m_newSeq_Name.find(seq);
            if (it != m_newSeq_Name.end()) {
                NS_LOG_DEBUG("Erase the tracing in retransmission. All data received from sequence number: " << seq);
                m_newSeq_Name.erase(seq);
            }

            // get mean average and serialize data
            /// !!!!!!Currently only don't compute average at aggregator, only add them together, perform average at consumer!!!!
            std::vector<uint8_t> newbuffer;
            serializeModelData(getMean(seqNum), newbuffer);

            // create data packet
            auto data = make_shared<Data>();

            std::string dest;

            dest = m_agg_newDataName[seq];
            std::string name_string = m_prefix.toUri() + "/" + dest + "/data" + "/" + seqNum;

            NS_LOG_DEBUG("New aggregated data's name: " << name_string);

            shared_ptr<Name> newName = make_shared<Name>(name_string);
            data->setName(*newName);
            data->setContent(make_shared< ::ndn::Buffer>(newbuffer.begin(), newbuffer.end()));
            data->setFreshnessPeriod(::ndn::time::milliseconds(m_freshness.GetMilliSeconds()));
            SignatureInfo signatureInfo(static_cast< ::ndn::tlv::SignatureTypeValue>(255));

            if (m_keyLocator.size() > 0){
                signatureInfo.setKeyLocator(m_keyLocator);
            }

            data->setSignatureInfo(signatureInfo);
            ::ndn::EncodingEstimator estimator;
            ::ndn::EncodingBuffer encoder(estimator.appendVarNumber(m_signature), 0);
            encoder.appendVarNumber(m_signature);
            data->setSignatureValue(encoder.getBuffer());

            // Send data packet
            data->wireEncode();
            m_transmittedDatas(data, this, m_face);
            m_appLink->onReceiveData(*data);
        } else{
            NS_LOG_DEBUG("Wait for others to aggregate.");
        }
    }else{
        NS_LOG_DEBUG("error when OnData()!");
    }
}

} // namespace ndn
} // namespace ns3


