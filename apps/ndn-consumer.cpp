/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/**
 * Copyright (c) 2011-2015  Regents of the University of California.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors and
 * contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

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

#include "utils/ndn-ns3-packet-tag.hpp"
#include "utils/ndn-rtt-mean-deviation.hpp"

#include <ndn-cxx/lp/tags.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/ref.hpp>

NS_LOG_COMPONENT_DEFINE("ndn.Consumer");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(Consumer);

TypeId
Consumer::GetTypeId(void)
{
  static TypeId tid =
    TypeId("ns3::ndn::Consumer")
      .SetGroupName("Ndn")
      .SetParent<App>()
      .AddAttribute("StartSeq", "Initial sequence number", IntegerValue(0),
                    MakeIntegerAccessor(&Consumer::m_seq), MakeIntegerChecker<int32_t>())
      .AddAttribute("Prefix", "Name of the Interest", StringValue("/"),
                    MakeNameAccessor(&Consumer::m_interestName), MakeNameChecker())
      .AddAttribute("LifeTime", "LifeTime for interest packet", StringValue("2s"),
                    MakeTimeAccessor(&Consumer::m_interestLifeTime), MakeTimeChecker())

      .AddAttribute("RetxTimer",
                    "Timeout defining how frequent retransmission timeouts should be checked",
                    StringValue("50ms"),
                    MakeTimeAccessor(&Consumer::GetRetxTimer, &Consumer::SetRetxTimer),
                    MakeTimeChecker())

      .AddTraceSource("LastRetransmittedInterestDataDelay",
                      "Delay between last retransmitted Interest and received Data",
                      MakeTraceSourceAccessor(&Consumer::m_lastRetransmittedInterestDataDelay),
                      "ns3::ndn::Consumer::LastRetransmittedInterestDataDelayCallback")

      .AddTraceSource("FirstInterestDataDelay",
                      "Delay between first transmitted Interest and received Data",
                      MakeTraceSourceAccessor(&Consumer::m_firstInterestDataDelay),
                      "ns3::ndn::Consumer::FirstInterestDataDelayCallback");

  return tid;
}

Consumer::Consumer()
    : m_rand(CreateObject<UniformRandomVariable>())       // Used to set Nonce
    , globalRound(0)
    , globalSeq(0)
    , SRTT(0)
    , RTTVAR(0)
    , roundRTT(0)
    , total_response_time(0)
    , round(0)
    , totalAggregateTime(0)
    , iteration(0)
    , m_seq(0)                                            // Used to set sequence num to 0 first, and keep a count of interest and manage interest names
{
    m_rtt = CreateObject<RttMeanDeviation>();
}


std::vector<std::string>
Consumer::getProducers()
{
    std::ifstream file(filename);
    std::vector<std::string> proNodes;
    std::string line;
    bool inRouterSection = false;

    if (!file) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return proNodes; // Return empty if file cannot be opened
    }

    while (std::getline(file, line)) {
        // Trim whitespace and check if the line is a section header
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v")); // Trim leading whitespace
        if (line == "router") {
            inRouterSection = true;
        } else if (line == "link") {
            break; // Exit if we reach the "link" section
        } else if (inRouterSection) {
            // Check if the line starts with "pro"
            if (line.substr(0, 3) == "pro") {
                proNodes.push_back(line);
            }
        }
    }
    file.close();
    return proNodes;
}


void
Consumer::AggregateTimeSum(int64_t aggregate_time){
    totalAggregateTime += aggregate_time;
    NS_LOG_DEBUG("totalAggregateTime is: " << totalAggregateTime);
    ++iteration;
}

int64_t
Consumer::GetAggregateTimeAverage() {
    if (iteration == 0) {
        NS_LOG_DEBUG("Error when calculating aggregate time!");
        return 0;
    }
    return totalAggregateTime / iteration;
}

void
Consumer::ResponseTimeSum(int64_t response_time)
{
    total_response_time += response_time;
    ++round;
}



void
Consumer::SetRetxTimer(Time retxTimer)
{
    m_retxTimer = retxTimer;
    if (m_retxEvent.IsRunning()) {
        Simulator::Remove(m_retxEvent); // slower, but better for memory
    }

    // Schedule new timeout
    m_timeoutThreshold = 6 * m_retxTimer;
    NS_LOG_DEBUG("Next interval to check timeout is: " << m_retxTimer.GetMilliSeconds() << " ms");
    m_retxEvent = Simulator::Schedule(m_retxTimer, &Consumer::CheckRetxTimeout, this);
}

Time
Consumer::GetRetxTimer() const
{
  return m_retxTimer;
}


void
Consumer::CheckRetxTimeout()
{
    Time now = Simulator::Now();

    NS_LOG_DEBUG("Check timeout after: " << m_retxTimer.GetMilliSeconds() << " ms");
    NS_LOG_DEBUG("Current timeout threshold is: " << m_timeoutThreshold.GetMilliSeconds() << " ms");

    //NS_LOG_DEBUG("Start timeout check:");
    for (auto it = m_timeoutCheck.begin(); it != m_timeoutCheck.end();){
        if (now - it->second > m_timeoutThreshold) {
            //NS_LOG_DEBUG("Packet " << it->first << " reach threshold, timeout triggered.");
            std::string name = it->first;
            it = m_timeoutCheck.erase(it);
                OnTimeout(name);
        } else {
            //NS_LOG_DEBUG("No timeout for packet " << it->first);
            ++it;
        }
    }
    m_retxEvent = Simulator::Schedule(m_retxTimer, &Consumer::CheckRetxTimeout, this);
}

/**
 * Compute new RTO based on response time of recent packets
 * @param resTime
 * @return New RTO
 */
Time
Consumer::RTOMeasurement(int64_t resTime)
{
    if (roundRTT == 0) {
        RTTVAR = resTime / 2;
        SRTT = resTime;
    } else {
        RTTVAR = 0.75 * RTTVAR + 0.25 * std::abs(SRTT - resTime); // RTTVAR = (1 - b) * RTTVAR + b * |SRTT - RTTsample|, where b = 0.25
        SRTT = 0.875 * SRTT + 0.125 * resTime; // SRTT = (1 - a) * SRTT + a * RTTsample, where a = 0.125
    }
    roundRTT++;
    int64_t RTO = SRTT + 4 * RTTVAR; // RTO = SRTT + K * RTTVAR, where K = 4

    return MilliSeconds(RTO);
}


// Application Methods
void
Consumer::StartApplication() // Called at time specified by Start
{
    //NS_LOG_FUNCTION_NOARGS();

    // Clear the log file
    std::ofstream file1(RTO_recorder, std::ios::out);
    if (!file1.is_open()) {
        std::cerr << "Failed to open the file: " << RTO_recorder << std::endl;
    }
    file1.close();

    std::ofstream file2(responseTime_recorder, std::ios::out);
    if (!file2.is_open()) {
        std::cerr << "Failed to open the file: " << responseTime_recorder << std::endl;
    }
    file2.close();

    Simulator::Schedule(MilliSeconds(5), &Consumer::RTORecorder, this);

    // do base stuff
    App::StartApplication();
    proNames = getProducers();
    for (const auto& item : proNames)
        NS_LOG_INFO("proNames: " << item);
    if (proNames.empty())
        NS_LOG_DEBUG("Error happened when getting producers!!!");

    ScheduleNextPacket();
}

void
Consumer::StopApplication() // Called at time specified by Stop
{
  NS_LOG_FUNCTION_NOARGS();

  // cancel periodic packet generation
  Simulator::Cancel(m_sendEvent);

  // cleanup base stuff
  App::StopApplication();
}

void
Consumer::SendPacket()
{
  if (!m_active)
    return;

  // Aggregation time start
  if (globalRound == 0) {
      aggregateStartTime[globalSeq] = ns3::Simulator::Now();
  }

  if (globalRound != proNames.size()) {
      std::string nameString = "/" + proNames[globalRound];
      map_agg[globalSeq].push_back(proNames[globalRound]);
      shared_ptr<Name> name = make_shared<Name>(nameString);
      name->appendSequenceNumber(globalSeq);

      globalRound++;
      SendInterest(name);
  } else {
      globalRound = 0;
      globalSeq = ++m_seq;
  }
  ScheduleNextPacket();

}


void
Consumer::SendInterest(shared_ptr<Name> newName) {

    std::string nameWithSeq = newName->toUri();
    // Trace timeout
    m_timeoutCheck[nameWithSeq] = ns3::Simulator::Now();

    // Start response time
    if (currentTime.find(nameWithSeq) == currentTime.end())
        currentTime[nameWithSeq] = ns3::Simulator::Now();

    shared_ptr<Interest> interest = make_shared<Interest>();
    interest->setNonce(m_rand->GetValue(0, std::numeric_limits<uint32_t>::max()));
    interest->setName(*newName);
    interest->setCanBePrefix(false);
    time::milliseconds interestLifeTime(m_interestLifeTime.GetMilliSeconds());
    interest->setInterestLifetime(interestLifeTime);
    NS_LOG_INFO("Sending interest >>>> " << nameWithSeq);
    m_transmittedInterests(interest, this, m_face);
    m_appLink->onReceiveInterest(*interest);
}


///////////////////////////////////////////////////
//          Process incoming packets             //
///////////////////////////////////////////////////

void
Consumer::OnData(shared_ptr<const Data> data)
{
  if (!m_active)
    return;

  App::OnData(data); // tracing inside

  //NS_LOG_FUNCTION(this << data);

  NS_LOG_INFO ("Received content object: " << boost::cref(*data));
  NS_LOG_INFO("The incoming data packet size is: " << data->wireEncode().size());

  std::string dataName = data->getName().toUri();
  uint32_t seq = data->getName().at(-1).toSequenceNumber();

  m_timeoutCheck.erase(dataName);
  NS_LOG_DEBUG("Erase packet: " << dataName << " from timeout list!");

  /// Customized
  std::string name = data->getName().get(0).toUri();
  auto data_map = map_agg.find(seq);
  if (data_map != map_agg.end()) {
      auto& vec = data_map->second;
      auto vecIt = std::find(vec.begin(), vec.end(), name);

      if (vecIt != vec.end()) {
          vec.erase(vecIt);
      } else {
          NS_LOG_INFO("Data name doesn't exist!");
          return;
      }

      // RTT computation
      if (currentTime.find(dataName) != currentTime.end()) {
          responseTime[dataName] = ns3::Simulator::Now() - currentTime[dataName];
          ResponseTimeSum(responseTime[dataName].GetMilliSeconds());
          currentTime.erase(dataName);
          NS_LOG_INFO("Consumer's response time of sequence " << dataName << " is: " << responseTime[dataName].GetMilliSeconds() << " ms");
      }

      // Record response time
      responseTimeRecorder(responseTime[dataName]);

      // Reset RetxTimer and timeout interval
      RTO_Timer = RTOMeasurement(responseTime[dataName].GetMilliSeconds());
      NS_LOG_DEBUG("responseTime for name : " << dataName << " is: " << responseTime[dataName].GetMilliSeconds() << " ms");
      NS_LOG_DEBUG("RTT measurement: " << RTO_Timer.GetMilliSeconds() << " ms");
      m_timeoutThreshold = RTO_Timer;

      if (vec.empty()) {
          NS_LOG_DEBUG("Aggregation iteration " << std::to_string(seq) << " finished!");

          if (aggregateStartTime.find(seq) != aggregateStartTime.end()) {
              aggregateTime[seq] = ns3::Simulator::Now() - aggregateStartTime[seq];
              AggregateTimeSum(aggregateTime[seq].GetMilliSeconds());
              aggregateStartTime.erase(seq);

          } else {
              NS_LOG_DEBUG("Error when calculating aggregation time!");
          }
      }

      /// Stop simulation
      if (iteration == 100) {
          NS_LOG_DEBUG("Reach 100 iterations, stop!");
          ns3::Simulator::Stop();
          NS_LOG_INFO("The average aggregation time of Consumer in " << iteration << " iteration is: " << GetAggregateTimeAverage() << " ms");
          return;
      }

  } else {
      NS_LOG_DEBUG("Packet " << dataName << " is not in the list! Error!!!!!");
  }

}

void
Consumer::OnNack(shared_ptr<const lp::Nack> nack)
{
  /// tracing inside
  App::OnNack(nack);

  NS_LOG_INFO("NACK received for: " << nack->getInterest().getName()
              << ", reason: " << nack->getReason());
}

void
Consumer::OnTimeout(std::string nameString)
{
    NS_LOG_DEBUG(nameString << " has timeout.");
    shared_ptr<Name> name = make_shared<Name>(nameString);
    SendInterest(name);
}


/**
 * Record the RTT every 5 ms, store them in a file
 */
void
Consumer::RTORecorder()
{
    // Open the file using fstream in append mode
    std::ofstream file(RTO_recorder, std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << RTO_recorder << std::endl;
        return;
    }

    // Write the response_time to the file, followed by a newline
    file << ns3::Simulator::Now().GetMilliSeconds() << " " << RTO_Timer.GetMilliSeconds() << std::endl;

    // Close the file
    file.close();
    Simulator::Schedule(MilliSeconds(5), &Consumer::RTORecorder, this);
}



/**
 * Record the response time for each returned packet, store them in a file
 * @param responseTime
 */
void
Consumer::responseTimeRecorder(Time responseTime) {
    // Open the file using fstream in append mode
    std::ofstream file(responseTime_recorder, std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << responseTime_recorder << std::endl;
        return;
    }

    // Write the response_time to the file, followed by a newline
    file << ns3::Simulator::Now().GetMilliSeconds() << " " << responseTime.GetMilliSeconds() << std::endl;

    // Close the file
    file.close();
}


} // namespace ndn
} // namespace ns3
